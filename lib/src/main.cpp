#include "gdal.h"
#include "GDAL_OPENCV_IO.h"
#include "spatial_datasets/patchified_w_raster_dataset.h"
#include "spatial_datasets/patchified_r_raster_dataset.h"
#include "io_raster.h"
#include "io_shapefile.h"
//#include "fusion_ops/epipolar_disparities_combination.h"
//#include "sample_stereo_matching.h"
#include "contours_stitching_evaluation.h"
#include "sample_zonal_stats.h"
#include "probability_to_proximity.h"

using namespace LxGeo::LxSpatialOps;
using namespace LxGeo::GeometryFactoryShared;

int main(int argc, char* argv[])
{
	clock_t t_begin = clock();
	GDALAllRegister();

	//std::string ref_raster_map = "C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/pseudo_dsm_inv.tif";
	//std::string in_polygons = "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/preds/build_poly.shp";
	//std::string out_polygons = "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/preds/build_poly.shp";
	//ZonalStatsPolygon dd = ZonalStatsPolygon()

	clock_t t_end = clock();
	std::cout << "** Elapsed time : " << double(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}
