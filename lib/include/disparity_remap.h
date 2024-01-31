#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"

namespace LxGeo
{
	namespace LxSpatialOps
	{

		class DisparityRemapper : public BaseOperation {

		public:
			DisparityRemapper(
				const std::string& left_disparity_path,
				const std::string& out_left_disparity_path
			) {
				add_raster_input_dataset(left_disparity_path, "in_ld");
				add_raster_output_dataset(out_left_disparity_path, { "in_ld" }, ExtentsCombinationStrategy::ext_intersection, "out_ld");
				init_cpd(OperationDiveStrategy::zoom, 3000.0, 500.0);
			};

			ViewPair op(ViewPair& in_view_pair) override {

				int16_t TEMP_NO_DATA = -32768;

				ViewPair out_view;
				const cv::Mat& in_ld = in_view_pair.raster_views["in_ld"].image;

				// Init output maps
				out_view.raster_views["out_ld"] = GeoImage<cv::Mat>((cv::Mat)(cv::Mat::ones(in_ld.size(), CV_16SC1) * TEMP_NO_DATA), in_view_pair.raster_views["in_ld"].geotransform);
				out_view.raster_views["out_ld"].set_nodata(in_view_pair.raster_views["in_ld"].no_data);
				cv::Mat& out_left_image = out_view.raster_views["out_ld"].image;

				cv::Mat map_x(in_ld.size(), CV_32FC1);
				cv::Mat map_y(in_ld.size(), CV_32FC1);

				for (int i = 0; i < map_x.rows; i++)
				{
					for (int j = 0; j < map_x.cols; j++)
					{
						short vv = in_ld.at<short>(i, j);
						if (vv == TEMP_NO_DATA)
							continue;
						map_x.at<float>(i, j) =  j + float(vv) /10;
						map_y.at<float>(i, j) = (float)i;

					}
				}

				remap(in_ld, out_left_image, map_x, map_y, cv::INTER_NEAREST);

				return out_view;

			}
		};

	}
}