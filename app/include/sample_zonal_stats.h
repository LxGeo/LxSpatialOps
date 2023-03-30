#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"
#include "stitching/vector_on_raster_stitcher.h"
#include "numcpp/stats.h"

namespace LxGeo
{
	namespace LxRasterOps
	{

		class ZonalStatsPolygon : public BaseOperation {

		public:
			ZonalStatsPolygon(
				const std::string& reference_raster_path,
				const std::string& polygons_map_path,
				const std::string& out_polygons_map_path
			) {

				add_raster_input_dataset(reference_raster_path, "ref_raster");
				add_vector_input_dataset(polygons_map_path, "polys");
				add_vector_output_dataset(out_polygons_map_path, { "polys" }, ext_intersection, WriteMode::create, "out_polys");
				init_cpd(OperationDiveStrategy::zoom, 1000.0, 200.0);

				auto ref_profile = input_raster_defmaps["ref_raster"].second;
				null_value = (ref_profile.no_data.has_value()) ? ref_profile.no_data.value() : FLT_MAX;
			}

			float null_value;

			ViewPair op(ViewPair& in_view_pair) override {

				ViewPair out_view;

				out_view.vector_views["out_polys"] = boost::get<GeoVector<Boost_Polygon_2>>(in_view_pair.vector_views["polys"]);
				GeoVector<Boost_Polygon_2>& gvec = boost::get<GeoVector<Boost_Polygon_2>>(out_view.vector_views["out_polys"]);

				RasterPixelsStitcher rps = RasterPixelsStitcher(in_view_pair.raster_views["ref_raster"]);

				bool valid_geom_assigned = false;
				for (auto& gwa : gvec.geometries_container) {
					auto stitched_pixels = rps.readPolygonPixels<float>(gwa.get_definition(), RasterPixelsStitcherStartegy::filled_polygon);
					auto stats = numcpp::DetailedStats<float>(stitched_pixels, null_value, 0.0);
					if (!stats.empty()) {
						valid_geom_assigned = true;
						gwa.set_int_attribute("count", stats.count());
						gwa.set_int_attribute("count_null", stats.count_null());
						gwa.set_double_attribute("null_r", double(stats.count_null()) / stats.count());
						gwa.set_double_attribute("max", stats.max());
						gwa.set_double_attribute("min", stats.min());
						gwa.set_double_attribute("mean", stats.mean());
						gwa.set_double_attribute("variance", stats.variance());
						gwa.set_double_attribute("p90", stats.percentile(90));
					}
					else if (!valid_geom_assigned)
						out_view.valid_geometries_indices["out_polys"].idx += 1;
				}
				return out_view;

			}

		};

	}
}