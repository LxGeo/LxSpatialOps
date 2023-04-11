#include "gdal.h"
#include "parameters.h"
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
#include "multistage_alignement.h"

using namespace LxGeo::LxRasterOps;
using namespace LxGeo::GeometryFactoryShared;

int main(int argc, char* argv[])
{
	clock_t t_begin = clock();
	GDALAllRegister();
	
	// Runs process
	std::unordered_map<std::string, std::string> view1_paths_map = {
		{"ortho", "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/ade27182-11da-483d-8fef-fb1a76b00568_z19.tif"},
		{"proba", "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/preds/build_probas.tif"},
		{"dsm", "C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/pseudo_dsm_inv.tif"},
		{"vector", "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/Brazil_Vila_Velha_A_Neo.shp"}
	};
	std::unordered_map<std::string, std::string> view2_paths_map = {
		{"ortho", "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_B_Neo/b25caa5a-190b-4fa9-957e-43816cc462c2_z19.tif"},
		{"proba", "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_B_Neo/preds/build_probas.tif"},
		{"dsm", "C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/pseudo_dsm.tif"},
		{"vector", "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_B_Neo/Brazil_Vila_Velha_B_Neo.shp"}
	};

	std::string couple_file_path = "C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/couple.json";

	std::string out_vector1 = "C:/Users/geoimage/Pictures/v1.shp";
	std::string out_vector2 = "C:/Users/geoimage/Pictures/v2.shp";


	MultiStageAlignment dd = MultiStageAlignment(view1_paths_map, view2_paths_map, couple_file_path, out_vector1, out_vector2);
	dd.run_sequential();

	//std::string ref_raster_map = "C:/DATA_SANDBOX/Alignment_Project/sgbm_pipeline/perfect_gt_brazil/pseudo_dsm_inv.tif";
	//std::string in_polygons = "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/preds/build_poly.shp";
	//std::string out_polygons = "C:/DATA_SANDBOX/Alignment_Project/PerfectGT/Brazil_Vila_Velha_A_Neo/preds/build_poly.shp";
	//ZonalStatsPolygon dd = ZonalStatsPolygon()

	clock_t t_end = clock();
	std::cout << "** Elapsed time : " << double(t_end - t_begin) / CLOCKS_PER_SEC << " s." << std::endl;

	return 0;
}
