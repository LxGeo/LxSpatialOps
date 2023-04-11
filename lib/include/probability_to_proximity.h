#pragma once
#include "base_operation.h"
#include "gdal_algs_wrap/gdal_proximity.h"
#include "defs.h"

namespace LxGeo
{
	namespace LxSpatialOps
	{

		class Proba2Prox : public BaseOperation {

		public:
			Proba2Prox(
				const std::string& input_proba,
				const std::string& output_proximity,
				const std::string& height_field_name_v2 = "h_median"
			) {

				add_raster_input_dataset(input_proba, "input_proba");
				add_raster_output_dataset(output_proximity, { "input_proba" }, ExtentsCombinationStrategy::ext_intersection, "output_proximity");
				init_cpd(OperationDiveStrategy::zoom, 1000.0, 200.0);

			}

			ViewPair op(ViewPair& in_view_pair) override {
				ViewPair out_view;

				const GeoImage<cv::Mat>& proba = in_view_pair.raster_views["input_proba"];
				const std::function<cv::Mat (cv::Mat)> binarizer_fn = [](const cv::Mat& in_img) {
					
					cv::Mat out_img = cv::Mat::zeros(in_img.size(), CV_8UC1)* in_img.channels();
					std::vector<cv::Mat> planes;
					cv::split(in_img, planes);

					for (int b_idx = 0; b_idx < planes.size(); b_idx++) {
						cv::Mat b_is_larger = cv::Mat::ones(in_img.size(), CV_8UC1);
						for (int b_idx2 = 0; b_idx2 < planes.size(); b_idx2++) {
							cv::Mat binaryImage;
							cv::compare(planes[b_idx], planes[b_idx2], binaryImage, cv::CMP_GE);
							cv::bitwise_and(b_is_larger, binaryImage, b_is_larger);
						}
						out_img += b_is_larger * (b_idx);
					}
					return out_img;

					//cv::reduceArgMax(in_img, out_img, 2);
					
				};
				VirtualGeoImage<cv::Mat> binary_contour = VirtualGeoImage<cv::Mat>(proba, binarizer_fn);
				GeoImage<cv::Mat> proximity_gimg = proximity_raster(binary_contour, 1, {2});
				out_view.raster_views["output_proximity"] = proximity_gimg;
				return out_view;
			}


		};
	}
} 