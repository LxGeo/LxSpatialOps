#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"

namespace LxGeo
{
	namespace LxRasterOps
	{

		class StereoMatching : public BaseOperation {

		public:
			StereoMatching(
				const std::string& left_image_path,
				const std::string& right_image_path,
				const std::string& out_disparity_path
			) {
				add_raster_input_dataset(left_image_path, "in_li");
				add_raster_input_dataset(right_image_path, "in_ri");
				add_raster_output_dataset(out_disparity_path, { "in_li" }, ExtentsCombinationStrategy::ext_intersection, "out_d");
				init_cpd(OperationDiveStrategy::zoom, 1000.0, 200.0);
			};

			ViewPair op(ViewPair& in_view_pair) override {

				ViewPair out_view;
				const cv::Mat& in_li = in_view_pair.raster_views["in_li"].image;
				const cv::Mat& in_ri = in_view_pair.raster_views["in_ri"].image;

				std::vector<cv::Mat> bands_left; cv::split(in_li, bands_left);
				std::vector<cv::Mat> bands_right; cv::split(in_ri, bands_right);
				bands_left[0].convertTo(bands_left[0], CV_8U);
				bands_right[0].convertTo(bands_right[0], CV_8U);

				cv::cuda::GpuMat in_li_g, in_ri_g, out_d_g, filtered_d_g;
				in_li_g.upload(bands_left[0]); 
				in_ri_g.upload(bands_right[0]);

				cv::Mat out_d;
				auto sgbm = cv::cuda::createStereoSGM();
				auto disparity_bilateral_filter = cv::cuda::createDisparityBilateralFilter(sgbm->getNumDisparities());

				sgbm->compute(in_li_g, in_ri_g, out_d_g);
				//out_d_g.download(out_d);

				disparity_bilateral_filter->apply(out_d_g, in_li_g, filtered_d_g);
				filtered_d_g.download(out_d);
				// Init output maps
				out_view.raster_views["out_d"] = GeoImage<cv::Mat>(out_d, in_view_pair.raster_views["in_li"].geotransform);
				
				return out_view;

			}
		};

	}
}