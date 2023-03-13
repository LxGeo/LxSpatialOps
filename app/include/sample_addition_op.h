#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"

namespace LxGeo
{
	namespace LxRasterOps
	{

		class MapAddition : public BaseOperation {

		public:
			MapAddition(
				const std::string& im1,
				const std::string& im2,
				const std::string& out_im
			) {

				add_raster_input_dataset(im1, "im1");
				add_raster_input_dataset(im2, "im2");
				add_raster_output_dataset(out_im, {"im1"}, ExtentsCombinationStrategy::ext_intersection, "out_im");
				init_cpd(OperationDiveStrategy::same);
			};

			ViewPair op(ViewPair& in_view_pair) override {
				ViewPair out_view;
				auto& r_views = in_view_pair.raster_views;
				const GeoImage<cv::Mat>& gimg1 = r_views["im1"];
				const GeoImage<cv::Mat>& gimg2 = r_views["im2"];
				cv::Mat out_arr = gimg1.image + gimg1.image;
				GeoImage<cv::Mat> out_img(out_arr, gimg1.geotransform);

				out_view.raster_views["out_im"] = out_img;
				return out_view;
			}
		};

	}
}