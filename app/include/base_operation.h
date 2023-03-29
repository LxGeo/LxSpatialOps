#pragma once
#include "defs.h"
#include "spatial_datasets/patchified_dataset.h"
#include "lightweight\raster_profile.h"
#include "lightweight\geoimage.h"
#include "spatial_datasets/patchified_w_raster_dataset.h"

#include "lightweight\geovector.h"
#include "lightweight\vector_profile.h"
#include "spatial_datasets/patchified_w_vector_dataset.h"

namespace LxGeo
{
	namespace LxRasterOps 
	{

		using namespace IO_DATA;

		enum ExtentsCombinationStrategy {
			ext_intersection = 0,
			ext_union = 1 << 0
		};
		Boost_Box_2 combine_extents(const std::list<Boost_Box_2>& input_extents_list, ExtentsCombinationStrategy ecs) {

			assert(input_extents_list.size() > 0);
			auto c_box = input_extents_list.begin();
			double o_minX= c_box->min_corner().get<0>(), o_maxX= c_box->max_corner().get<0>(), o_minY= c_box->min_corner().get<1>(), o_maxY= c_box->max_corner().get<1>();
			std::advance(c_box, 1);

			for (; c_box != input_extents_list.end(); c_box++) {
				if (ecs == ExtentsCombinationStrategy::ext_intersection) {
					o_minX = std::max(o_minX, c_box->min_corner().get<0>());
					o_maxX = std::min(o_maxX, c_box->max_corner().get<0>());
					o_minY = std::max(o_minY, c_box->min_corner().get<1>());
					o_maxY = std::min(o_maxY, c_box->max_corner().get<1>());
				}
				else if (ecs == ExtentsCombinationStrategy::ext_union) {
					o_minX = std::min(o_minX, c_box->min_corner().get<0>());
					o_maxX = std::max(o_maxX, c_box->max_corner().get<0>());
					o_minY = std::min(o_minY, c_box->min_corner().get<1>());
					o_maxY = std::max(o_maxY, c_box->max_corner().get<1>());
				}
				else {
					throw std::exception("Not emplemented yet");
				}
			}
			return Boost_Box_2({ o_minX, o_minY }, { o_maxX, o_maxY });
		}

		// Strategy used to generate continuous patches for input and output datasets.
		enum OperationDiveStrategy {
			same = 0,
			zoom = 1 <<0
		};

		typedef boost::variant< GeoVector<Boost_Polygon_2>, GeoVector<Boost_LineString_2>, GeoVector<Boost_Point_2> > GeoVectorVariant;
		struct ViewPair {
			std::map<std::string, GeoImage<cv::Mat>> raster_views;
			std::map<std::string, GeoVectorVariant> vector_views;
		};

		/**
		Base class for spatial dataset operation. Defines requierement for operations.
		*/
		class BaseOperation{
		
		public:
			void add_raster_input_dataset(std::string dataset_path, std::string dataset_id = "") {
				assert(input_raster_defmaps.find(dataset_id) == input_raster_defmaps.end() && "Dataset ID already defined!");
				auto dst_ptr = load_gdal_dataset_shared_ptr(dataset_path);
				auto c_dst_profile = IO_DATA::RProfile::from_gdal_dataset(dst_ptr);
				dataset_id = (dataset_id.empty()) ? dataset_path : dataset_id;
				input_raster_defmaps[dataset_id] = { dataset_path, c_dst_profile };
			}
			void add_raster_output_dataset(
				std::string dataset_path,  std::set<std::string> respective_datasets_ids,
				ExtentsCombinationStrategy ecs= ExtentsCombinationStrategy::ext_intersection, std::string dataset_id = ""
			) {
				assert(respective_datasets_ids.size() > 0, "Respective datasets IDs is empty!");
				const IO_DATA::RProfile& first_profile = input_raster_defmaps[*respective_datasets_ids.begin()].second;

				std::list<Boost_Box_2> extents;
				for (const auto& key_value : input_raster_defmaps) {
					if (respective_datasets_ids.find(key_value.first) == respective_datasets_ids.end())
						continue;
					auto c_dataset_profile = key_value.second.second;
					assert(
						c_dataset_profile.compare(first_profile, RProfileCompFlags::gsd),
						fmt::format("Dataset don't share the same gds! Check rasters at: \n {}", key_value.second.first)
					);
					assert(
						c_dataset_profile.compare(first_profile, RProfileCompFlags::crs_wkt),
						fmt::format("Dataset don't share the same spatial reference system! Check rasters at: \n {}", key_value.second.first)
					);
					extents.push_back(c_dataset_profile.to_box_extents());

				}

				Boost_Box_2 output_extents = combine_extents(extents, ecs);
				RProfile output_profile = RProfile::from_extents(output_extents, first_profile.geotransform[1], first_profile.geotransform[5], 1, GDT_Byte, const_cast<char*>(first_profile.s_crs_wkt.c_str()));

				dataset_id = (dataset_id.empty()) ? dataset_path : dataset_id;
				output_raster_defmaps[dataset_id] = { dataset_path, output_profile };

			}

			void add_vector_input_dataset(std::string dataset_path, std::string dataset_id = "") {
				assert(input_vector_defmaps.find(dataset_id) == input_vector_defmaps.end() && "Dataset ID already defined!");
				auto dst_ptr = load_gdal_vector_dataset_shared_ptr(dataset_path);
				auto c_dst_profile = IO_DATA::VProfile::from_gdal_dataset(dst_ptr);
				dataset_id = (dataset_id.empty()) ? dataset_path : dataset_id;
				input_vector_defmaps[dataset_id] = { dataset_path, c_dst_profile };
			}

			/**
			* respective_datasets_ids is used to initialize output vector profile. If respective dataset id correspond to a raster id then spatial reference system is copied,
			* else if id corresponds to a vector id then the whole VProfile is copied including layer definition (layer definition is goining to be updated later at init_output_datasets)
			*/
			void add_vector_output_dataset(std::string dataset_path, std::set<std::string> respective_datasets_ids,
				ExtentsCombinationStrategy ecs = ExtentsCombinationStrategy::ext_intersection, WriteMode wm = WriteMode::create, std::string dataset_id = "") {
				dataset_id = (dataset_id.empty()) ? dataset_path : dataset_id;

				std::list<IO_DATA::VProfile> v_profiles;
				std::list<IO_DATA::RProfile> r_profiles;
				std::list<Boost_Box_2> extents;
				std::unordered_set<std::string> spatial_refs;

				assert(respective_datasets_ids.size() > 0 && "Respective datasets IDs is empty!");
				for (const auto& c_respective_id : respective_datasets_ids) {
					if (input_raster_defmaps.find(c_respective_id) != input_raster_defmaps.end()) {
						auto& c_raster_profile = input_raster_defmaps[c_respective_id].second;
						r_profiles.push_back(c_raster_profile);
						extents.push_back(c_raster_profile.to_box_extents());
						spatial_refs.insert(c_raster_profile.s_crs_wkt);
					}
					else if (input_vector_defmaps.find(c_respective_id) != input_vector_defmaps.end()) {
						auto& c_vector_profile = input_vector_defmaps[c_respective_id].second;
						v_profiles.push_back(c_vector_profile);
						std::transform(
							c_vector_profile.layers_def.begin(), c_vector_profile.layers_def.end(),
							std::back_inserter(extents), [](const auto& layers_kv)->Boost_Box_2 {return layers_kv.second.layer_extent; }
						);
						spatial_refs.insert(c_vector_profile.s_crs_wkt);
					}
					else {
						auto err_msg = fmt::format("Respective ID: {} is not found in vector and raster definition maps!", c_respective_id);
						throw std::exception(err_msg.c_str());
					}
				}
				
				assert(spatial_refs.size() <= 1 && fmt::format("Respective datasets of ouput vector dataset with ID: {} don't share the same spatial reference system!", dataset_id).c_str());				
				Boost_Box_2 output_extents = combine_extents(extents, ecs);


				VProfile out_profile; out_profile.s_crs_wkt = *spatial_refs.begin();
				LayerDef ld; ld.layer_extent = output_extents;
				out_profile.layers_def[""] = ld;
				output_vector_defmaps[dataset_id] = { dataset_path, out_profile, wm };
			}
			
			void init_cpd(OperationDiveStrategy ods, double spatial_patch_size=256.0, double spatial_buffer = 0.0) {				
				std::list<Boost_Box_2> outputs_extents;
				for (auto& c_elem : output_raster_defmaps)
					outputs_extents.push_back(c_elem.second.second.to_box_extents());
				for (auto& c_elem : output_vector_defmaps)
					for (auto& c_layer : std::get<1>(c_elem.second).layers_def)
						outputs_extents.push_back(c_layer.second.layer_extent);
				
				// Computing union area for all output datasets
				Boost_Polygon_2 output_area; bg::assign(output_area, *outputs_extents.begin());
				for (auto& c_extent : outputs_extents) {
					std::deque<Boost_Polygon_2> union_geoms;
					bg::correct(output_area); bg::correct(c_extent);
					bg::union_(output_area, c_extent, union_geoms);
					assert(union_geoms.size() == 1, "Union of extents generated multipart polygons!");
					output_area = union_geoms.front();
				}

				out_cpd = ContinuousPatchifiedDataset({ spatial_patch_size, 0.0, 0.0, output_area });
				if (out_cpd.length() < 1)
					throw std::exception("Misgenerated patch grids! Make sure that spatial patch size corresponds well to extents!");
				if (ods == OperationDiveStrategy::same) {
					in_cpd = ContinuousPatchifiedDataset(out_cpd);
				}
				else if (ods == OperationDiveStrategy::zoom) {
					assert(spatial_buffer > 0.0, "OperationDiveStrategy is zoom! spatial buffer should be heigher than zero!");
					in_cpd = ContinuousPatchifiedDataset({ spatial_patch_size + spatial_buffer, spatial_buffer, spatial_buffer / 2, output_area });
				}

				assert(out_cpd.length() == in_cpd.length(), "Error at generating in & out patch grids!");
			}

			/*ViewPair operator[](int offset) {
				ViewPair c_view_pair;
				const Boost_Box_2& c_view_box = in_cpd[offset];
				// Read raster inputs
				for (const auto& c_input_raster_kv : input_raster_defmaps) {
					auto c_raster_id = c_input_raster_kv.first;
					GeoImage<cv::Mat> c_geoimage = GeoImage<cv::Mat>::from_file(c_input_raster_kv.second.first, c_view_box);
					c_view_pair[c_raster_id] = c_geoimage;
				}
				return c_view_pair;
			}*/

			virtual ViewPair op(ViewPair& in_view_pair) = 0;

			void init_output_datasets_from_view_pair(const ViewPair& out_view_pair) {
				for (const auto& out_view_kv : out_view_pair.raster_views) {
					if (output_raster_defmaps.find(out_view_kv.first) == output_raster_defmaps.end()) {
						throw std::exception("Operation generated an undefined raster view ID!");
					}
					auto& respective_profile = output_raster_defmaps[out_view_kv.first].second;
					respective_profile.count = out_view_kv.second.image.channels();
					respective_profile.dtype = KGDAL2CV::opencv2gdal(out_view_kv.second.image.type());
					output_raster_datasets[out_view_kv.first] = WPRasterDataset(output_raster_defmaps[out_view_kv.first].first, respective_profile, WriteMode::overwrite);
				}

				for (const auto& out_view_kv : out_view_pair.vector_views) {
					if (output_vector_defmaps.find(out_view_kv.first) == output_vector_defmaps.end()) {
						throw std::exception("Operation generated an undefined vector view ID!");
					}
					auto& respective_profile = std::get<1>(output_vector_defmaps[out_view_kv.first]);
					const GeoVectorVariant& gvv = out_view_kv.second;
					LayerDef c_layer_def;
					if (const GeoVector<Boost_Polygon_2>* poly_geov = boost::get<GeoVector<Boost_Polygon_2>>(&gvv)) {
						if (poly_geov->length() == 0) {
							std::cout << fmt::format("Empty polygon geometry contatiner at GeoVector corresponding to output vector with ID: {}",out_view_kv.first)<<std::endl;
							std::cout << "Matbe try to use a different grid index to initialize output datasets" << std::endl;
							continue;
						}
						c_layer_def = LayerDef::from_geometry_wa(poly_geov->geometries_container[0]);
					}
					else if (const GeoVector<Boost_LineString_2>* line_geov = boost::get<GeoVector<Boost_LineString_2>>(&gvv)) {
						if (line_geov->length() == 0) {
							std::cout << fmt::format("Empty line geometry contatiner at GeoVector corresponding to output vector with ID: {}", out_view_kv.first) << std::endl;
							std::cout << "Matbe try to use a different grid index to initialize output datasets" << std::endl;
							continue;
						}
						c_layer_def = LayerDef::from_geometry_wa(line_geov->geometries_container[0]);
					}
					else if (const GeoVector<Boost_Point_2>* point_geov = boost::get<GeoVector<Boost_Point_2>>(&gvv)) {
						if (point_geov->length() == 0) {
							std::cout << fmt::format("Empty point geometry contatiner at GeoVector corresponding to output vector with ID: {}", out_view_kv.first) << std::endl;
							std::cout << "Matbe try to use a different grid index to initialize output datasets" << std::endl;
							continue;
						}
						c_layer_def = LayerDef::from_geometry_wa(point_geov->geometries_container[0]);
					}
					else
						throw std::exception("Cannot initialize output vector dataset! Boost variant variable is out the three supported geometry types.");
					respective_profile.layers_def[""] = c_layer_def;
					output_vector_datasets[out_view_kv.first] = WPVectorDataset(std::get<0>(output_vector_defmaps[out_view_kv.first]), respective_profile, std::get<2>(output_vector_defmaps[out_view_kv.first]));
				}
			}

			void run_unique_block(int c_grid_id, bool init_output_datasets=false) {
				ViewPair in_view_pair;
				const Boost_Box_2& in_view_box = in_cpd[c_grid_id];
				// Read raster inputs
				for (const auto& c_input_raster_kv : input_raster_defmaps) {
					auto c_raster_id = c_input_raster_kv.first;
					GeoImage<cv::Mat> c_geoimage = GeoImage<cv::Mat>::from_file(c_input_raster_kv.second.first, in_view_box);
					in_view_pair.raster_views[c_raster_id] = c_geoimage;
				}
				// Read vector inputs
				for (const auto& c_input_vector_kv : input_vector_defmaps) {
					auto c_vector_id = c_input_vector_kv.first;
					auto& first_layer = *c_input_vector_kv.second.second.layers_def.begin();
					if (first_layer.second.wkb_type == wkbPoint) {
						GeoVector<Boost_Point_2> c_geovector = GeoVector<Boost_Point_2>::from_file<Boost_Box_2>(c_input_vector_kv.second.first, first_layer.first, in_view_box);
						in_view_pair.vector_views[c_vector_id] = c_geovector;
					}
					if (first_layer.second.wkb_type == wkbLineString) {
						GeoVector<Boost_LineString_2> c_geovector = GeoVector<Boost_LineString_2>::from_file<Boost_Box_2>(c_input_vector_kv.second.first, first_layer.first, in_view_box);
						in_view_pair.vector_views[c_vector_id] = c_geovector;
					}
					if (first_layer.second.wkb_type == wkbPolygon) {
						GeoVector<Boost_Polygon_2> c_geovector = GeoVector<Boost_Polygon_2>::from_file<Boost_Box_2>(c_input_vector_kv.second.first, first_layer.first, in_view_box);
						in_view_pair.vector_views[c_vector_id] = c_geovector;
					}
					
				}

				ViewPair out_view_pair = op(in_view_pair);
				if (init_output_datasets)
					init_output_datasets_from_view_pair(out_view_pair);
				const Boost_Box_2& out_view_box = out_cpd[c_grid_id];
				// Write raster outputs
				for (const auto& c_out_dataset_kv : out_view_pair.raster_views) {
					std::string c_dataset_id = c_out_dataset_kv.first;
					GeoImage<cv::Mat> respective_geoimage = out_view_pair.raster_views[c_dataset_id].get_view_spatial(out_view_box);
					output_raster_datasets[c_dataset_id].write_geoimage(respective_geoimage);
				}
				for (const auto& c_out_dataset_kv : out_view_pair.vector_views) {
					std::string c_dataset_id = c_out_dataset_kv.first;

					OGRwkbGeometryType c_dataset_type = std::get<1>(output_vector_defmaps[c_dataset_id]).layers_def.begin()->second.wkb_type;
					if (c_dataset_type == wkbPoint) {
						const GeoVector<Boost_Point_2>& respective_geovector = boost::get<GeoVector<Boost_Point_2>>(out_view_pair.vector_views[c_dataset_id]);
						output_vector_datasets[c_dataset_id].write_geovector(respective_geovector.get_view_spatial(out_view_box));
					}
					if (c_dataset_type == wkbLineString) {
						const GeoVector<Boost_LineString_2>& respective_geovector = boost::get<GeoVector<Boost_LineString_2>>(out_view_pair.vector_views[c_dataset_id]);
						output_vector_datasets[c_dataset_id].write_geovector(respective_geovector.get_view_spatial(out_view_box));
					}
					if (c_dataset_type == wkbPolygon) {
						const GeoVector<Boost_Polygon_2>& respective_geovector = boost::get<GeoVector<Boost_Polygon_2>>(out_view_pair.vector_views[c_dataset_id]);
						output_vector_datasets[c_dataset_id].write_geovector(respective_geovector.get_view_spatial(out_view_box));
					}

				}

			}

			void run_sequential() {
				run_unique_block(0, true);
				for (int c_grid_id = 1; c_grid_id < out_cpd.length(); c_grid_id++) {
					run_unique_block(c_grid_id);
				}
			}

			/*void run_parallel() {
				run_unique_block(0, true);
				std::vector<std::thread> threads; threads.reserve(out_cpd.length());
				for (int c_grid_id = 1; c_grid_id < out_cpd.length(); c_grid_id++) {
					threads.push_back(std::thread(&BaseOperation::run_unique_block, this, c_grid_id, false));
				}
				for (auto& th : threads) {
					th.join();
				}
			}*/

		public:

			ContinuousPatchifiedDataset out_cpd, in_cpd;
			std::map<std::string, std::pair<std::string, IO_DATA::RProfile>> input_raster_defmaps;
			std::map<std::string, std::pair<std::string, IO_DATA::RProfile>> output_raster_defmaps;

			std::map<std::string, std::pair<std::string, IO_DATA::VProfile>> input_vector_defmaps;
			std::map<std::string, std::tuple<std::string, IO_DATA::VProfile, WriteMode>> output_vector_defmaps;

			std::map<std::string, WPRasterDataset> output_raster_datasets;
			std::map<std::string, WPVectorDataset> output_vector_datasets;

			//std::map<> input_vector_datasets;
			//std::map<> output_vector_datasets;

		};

	}
}