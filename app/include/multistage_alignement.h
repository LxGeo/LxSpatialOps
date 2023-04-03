#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "gdal_algs_wrap/gdal_proximity.h"
#include "numcpp/stats.h"
#include <nlohmann/json.hpp>

namespace LxGeo
{
	namespace LxRasterOps
	{

		class MultiStageAlignment : public BaseOperation {

		public:
			std::pair<double, double> xy_constants;
			
		public:
			MultiStageAlignment(
				const std::unordered_map<std::string, std::string>& view1_paths_map,
				const std::unordered_map<std::string, std::string>& view2_paths_map,
				const std::string& couple_file_path,
				const std::string& aligned_vector_1,
				const std::string& aligned_vector_2
			) {

				add_raster_input_dataset(view1_paths_map.at("ortho"), "ortho_1");
				add_raster_input_dataset(view1_paths_map.at("proba"), "proba_1");
				add_raster_input_dataset(view1_paths_map.at("dsm"), "dsm_1");
				add_vector_input_dataset(view1_paths_map.at("vector"), "vector_1");

				add_raster_input_dataset(view2_paths_map.at("ortho"), "ortho_2");
				add_raster_input_dataset(view2_paths_map.at("proba"), "proba_2");
				add_raster_input_dataset(view2_paths_map.at("dsm"), "dsm_2");
				add_vector_input_dataset(view2_paths_map.at("vector"), "vector_2");

				add_vector_output_dataset(aligned_vector_1, { "vector_1" }, ext_intersection, WriteMode::overwrite, "al_v1");
				add_vector_output_dataset(aligned_vector_2, { "vector_2" }, ext_intersection, WriteMode::overwrite, "al_v2");

				std::ifstream f(couple_file_path);
				nlohmann::json data = nlohmann::json::parse(f);
				std::vector<double> loaded_constants = data["xy_cst"];
				if (loaded_constants.size() != 2)
					throw std::runtime_error("Error reading xy_constants from couple file");
				xy_constants = std::make_pair(loaded_constants[0], loaded_constants[1]);

				init_cpd(OperationDiveStrategy::zoom, 1000.0, 200.0);

			}

			GeoImage<cv::Mat> proba2proximity(const GeoImage<cv::Mat>& proba_gimg, size_t band_idx = 2) {
				// functor used to apply argmax on proba map and return a binary map of the respective band
				const std::function<cv::Mat(cv::Mat)> binarizer_fn = [&band_idx](const cv::Mat& in_img) {
					std::vector<cv::Mat> planes;
					cv::split(in_img, planes);
					cv::Mat b_is_larger = cv::Mat::ones(in_img.size(), CV_8UC1);
					for (int b_idx = 0; b_idx < planes.size(); b_idx++) {
						cv::Mat binaryImage;
						cv::compare(planes[band_idx], planes[b_idx], binaryImage, cv::CMP_GE);
						cv::bitwise_and(b_is_larger, binaryImage, b_is_larger);						
					}
					return b_is_larger;
				};
				VirtualGeoImage<cv::Mat> binary_contour = VirtualGeoImage<cv::Mat>(proba_gimg, binarizer_fn);
				GeoImage<cv::Mat> proximity_gimg = proximity_raster(binary_contour, 1);
				return proximity_gimg;
			}

			void assign_heights(GeoVector<Boost_Polygon_2>& in_g_vector, GeoImage<cv::Mat>& height_raster, const std::string& height_column_name = "height") {
				float null_value = (height_raster.no_data.has_value()) ? height_raster.no_data.value() : FLT_MAX; //value_or
				RasterPixelsStitcher rps = RasterPixelsStitcher(height_raster);
				for (auto& gwa : in_g_vector.geometries_container) {
					auto stitched_pixels = rps.readPolygonPixels<float>(gwa.get_definition(), RasterPixelsStitcherStartegy::filled_polygon);
					auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
					if (!stats.empty()) {
						gwa.set_double_attribute("null_r", double(stats.count_null()) / (stats.count()+ stats.count_null()));
						gwa.set_double_attribute(height_column_name, stats.percentile(90));
						gwa.set_int_attribute("H_LABEL", 1);
					}
					else {
						gwa.set_double_attribute("null_r", 0);
						gwa.set_double_attribute(height_column_name, DBL_MAX);
						gwa.set_int_attribute("H_LABEL", 0);
					}
				}
			}

			void assign_confidence(GeoVector<Boost_Polygon_2>& in_g_vector, GeoImage<cv::Mat>& proximity_raster, const std::string& confidence_column_name = "confidence") {
				float null_value = (proximity_raster.no_data.has_value()) ? proximity_raster.no_data.value() : FLT_MAX; //value_or
				RasterPixelsStitcher rps = RasterPixelsStitcher(proximity_raster);
				for (auto& gwa : in_g_vector.geometries_container) {
					auto stitched_pixels = rps.readPolygonPixels<float>(gwa.get_definition(), RasterPixelsStitcherStartegy::contours);
					auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
					if (!stats.empty()) {
						gwa.set_double_attribute(confidence_column_name, stats.sum()/stats.count());
						gwa.set_int_attribute("C_LABEL", 1);
					}
					else {
						gwa.set_double_attribute(confidence_column_name, DBL_MAX);
						gwa.set_int_attribute("C_LABEL", 0);
					}
				}
			}

			void roof2roof(GeoVector<Boost_Polygon_2>& in_g_vector, const std::pair<double, double>& xy_constants, const std::string& height_column_name = "height") {
				for (auto& gwa : in_g_vector.geometries_container) {
					int H_LABEL = gwa.get_int_attribute("H_LABEL");
					if (H_LABEL == 0)
						continue;
					double h = gwa.get_double_attribute(height_column_name);
					auto translate_matrix = bg::strategy::transform::translate_transformer<double, 2, 2>(h * xy_constants.first, h *xy_constants.second);
					gwa.set_definition(translate_geometry(gwa.get_definition(), translate_matrix));
				}
			};

			ViewPair op(ViewPair& in_view_pair) override {
				ViewPair out_view;

				GeoVector<Boost_Polygon_2>& polygons1 = boost::get<GeoVector<Boost_Polygon_2>>(in_view_pair.vector_views["vector_1"]);
				GeoVector<Boost_Polygon_2>& polygons2 = boost::get<GeoVector<Boost_Polygon_2>>(in_view_pair.vector_views["vector_2"]);

				assign_heights(polygons1, in_view_pair.raster_views["dsm_1"], "height");
				assign_heights(polygons2, in_view_pair.raster_views["dsm_2"], "height");

				roof2roof(polygons1, { -xy_constants.first, -xy_constants.second }, "height");
				roof2roof(polygons2, xy_constants, "height");

				GeoImage<cv::Mat> contour_proximity_1 = proba2proximity(in_view_pair.raster_views["proba_1"]);
				GeoImage<cv::Mat> contour_proximity_2 = proba2proximity(in_view_pair.raster_views["proba_2"]);

				assign_confidence(polygons1, contour_proximity_2, "confidence");
				assign_confidence(polygons2, contour_proximity_1, "confidence");

				out_view.valid_geometries_indices["al_v1"].idx = 0;
				out_view.valid_geometries_indices["al_v2"].idx = 0;
				out_view.vector_views["al_v1"] = std::move(polygons1);
				out_view.vector_views["al_v2"] = std::move(polygons2);

				return out_view;
			}


		};
	}
}