#include "gdal.h"
#include "parameters.h"
#include "GDAL_OPENCV_IO.h"
#include "spatial_datasets/patchified_w_raster_dataset.h"
#include "spatial_datasets/patchified_r_raster_dataset.h"
#include "io_raster.h"
#include "io_shapefile.h"
#include "fusion_ops/epipolar_disparities_combination.h"

using namespace LxGeo::LxRasterOps;
using namespace LxGeo::GeometryFactoryShared;

int main(int argc, char* argv[])
{
	clock_t t_begin = clock();
	GDALAllRegister();
	/*
	params = new Parameters(argc, argv);
	if (!params->initialized()) {
		delete params;
		return 1;
	}*/

	// Runs process
	RPRasterDataset<cv::Mat> s_dst= RPRasterDataset<cv::Mat>("C:/DATA_SANDBOX/lxProximityAlign/Brazil_Vila/sample1/temp_dir/proximity.tif", 256, 0);
	//WPRasterDataset d_dst = WPRasterDataset("C:/DATA_SANDBOX/lxProximityAlign/Brazil_Vila/sample1/temp_dir/proximity__.tif", s_dst.get_rprofile());
	//GeoImage<cv::Mat> a = s_dst[0];
	//d_dst.write_geoimage(s_dst[0]);

	DisparityMapsCombiner dd = DisparityMapsCombiner(
		"C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/disp.tif",
		"C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/disp_inv.tif",
		"C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/fusion/disp.tif",
		"C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/fusion/disp_inv.tif",
		"C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/fusion/disp_diff.tif");
	
	dd.run_sequential();
	// Quits
	delete params;

	clock_t t_end = clock();
	std::cout << "** Elapsed time : " << double(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}
