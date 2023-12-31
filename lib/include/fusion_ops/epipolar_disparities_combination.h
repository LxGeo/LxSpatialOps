#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"

namespace LxGeo
{
	namespace LxSpatialOps
	{

		class DisparityMapsCombiner: public BaseOperation {

		public:
			DisparityMapsCombiner(
				const std::string& left_disparity_path,
				const std::string& right_disparity_path,
				const std::string& out_left_disparity_path,
				const std::string& out_right_disparity_path,
				const std::string& left_certainty_map_path
			) {

				add_raster_input_dataset(left_disparity_path, "in_ld");
				add_raster_input_dataset(right_disparity_path, "in_rd");
				add_raster_output_dataset(out_left_disparity_path, {"in_ld"}, ExtentsCombinationStrategy::ext_intersection, "out_ld");
				add_raster_output_dataset(out_right_disparity_path, { "in_rd" }, ExtentsCombinationStrategy::ext_intersection, "out_rd");
				add_raster_output_dataset(left_certainty_map_path, { "in_ld" }, ExtentsCombinationStrategy::ext_intersection, "out_lcm");

				init_cpd(OperationDiveStrategy::zoom, 3000.0, 500.0);

			};

			ViewPair op (ViewPair& in_view_pair) override {

				int16_t TEMP_NO_DATA = -32768;

				ViewPair out_view;
				const cv::Mat& in_ld = in_view_pair.raster_views["in_ld"].image;
				const cv::Mat& in_rd = in_view_pair.raster_views["in_rd"].image;

				// Init output maps
				out_view.raster_views["out_ld"] = GeoImage<cv::Mat>((cv::Mat)(cv::Mat::ones(in_ld.size(), CV_16SC1)* TEMP_NO_DATA), in_view_pair.raster_views["in_ld"].geotransform); //GeoImage<cv::Mat>(in_view_pair.raster_views["in_ld"]);
				out_view.raster_views["out_rd"] = GeoImage<cv::Mat>((cv::Mat)(cv::Mat::ones(in_rd.size(), CV_16SC1) * TEMP_NO_DATA), in_view_pair.raster_views["in_rd"].geotransform);//GeoImage<cv::Mat>(in_view_pair.raster_views["in_rd"]);
				out_view.raster_views["out_lcm"] = GeoImage<cv::Mat>((cv::Mat)cv::Mat::zeros(in_ld.size(), CV_16SC1), in_view_pair.raster_views["in_ld"].geotransform);

				/*double l_px, r_px, diff;
				int xpos_in_other_image;
				for (int x = 0; x < in_ld.cols; x++) {
					for (int y = 0; y < in_ld.rows; y++) {
						// fill right image
						l_px = in_ld.at<int16_t>(y, x);
						xpos_in_other_image = x - l_px/10;
						
						if (xpos_in_other_image >= 0 && xpos_in_other_image < in_ld.cols && l_px!=TEMP_NO_DATA) {
							if (out_view.raster_views["out_rd"].image.at<int16_t>(y, xpos_in_other_image) == TEMP_NO_DATA)
								out_view.raster_views["out_rd"].image.at<int16_t>(y, xpos_in_other_image) = -l_px;
							else
								out_view.raster_views["out_lcm"].image.at<int16_t>(y,x) = abs(l_px + out_view.raster_views["out_rd"].image.at<int16_t>(y, xpos_in_other_image));
						}

						r_px = in_rd.at<int16_t>(y,x);
						xpos_in_other_image = x - r_px/10;
						if (xpos_in_other_image >= 0 && xpos_in_other_image < in_ld.cols && r_px!=TEMP_NO_DATA) {
							if (out_view.raster_views["out_ld"].image.at<int16_t>(y, xpos_in_other_image) == TEMP_NO_DATA)
								out_view.raster_views["out_ld"].image.at<int16_t>(y, xpos_in_other_image) = -r_px;
						}

					}
				}*/

				// Better only for release mode
				in_ld.forEach<int16_t>([&in_rd, &out_view, &TEMP_NO_DATA](int16_t& l_px, const int position[])->void {
					cv::Mat& out_rigth_image = out_view.raster_views["out_rd"].image;
					cv::Mat& out_left_image = out_view.raster_views["out_ld"].image;
					cv::Mat& out_left_cm = out_view.raster_views["out_lcm"].image;

					int xpos_in_other_image = position[1] - l_px / 10;
					if (xpos_in_other_image >= 0 && xpos_in_other_image < out_left_image.cols && l_px != TEMP_NO_DATA) {
						if (out_rigth_image.at<int16_t>(position[0], xpos_in_other_image) == TEMP_NO_DATA)
							out_rigth_image.at<int16_t>(position[0], xpos_in_other_image) = -l_px;
						else
							out_left_cm.at<int16_t>(position[0], position[1]) = abs(l_px + out_rigth_image.at<int16_t>(position[0], xpos_in_other_image));
					}

					int16_t r_px = in_rd.at<int16_t>(position[0], position[1]);
					xpos_in_other_image = position[1] - r_px / 10;
					if (xpos_in_other_image >= 0 && xpos_in_other_image < out_left_image.cols && r_px != TEMP_NO_DATA) {
						if (out_left_image.at<int16_t>(position[0], xpos_in_other_image) == TEMP_NO_DATA)
							out_left_image.at<int16_t>(position[0], xpos_in_other_image) = -r_px;
					}


					});


				return out_view;

			}
		};

	}
}