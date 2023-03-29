#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"
#include "stitching/vector_on_raster_stitcher.h"

namespace LxGeo
{
	namespace LxRasterOps
	{

		class ContoursStitchingEval : public BaseOperation {

		public:
			ContoursStitchingEval(
				const std::string& proximity_map_path,
				const std::string& polygons_map_path
			) {

				add_raster_input_dataset(proximity_map_path, "proximity");
				add_vector_input_dataset(polygons_map_path, "polys");
				add_vector_output_dataset(polygons_map_path, { "polys" }, ext_intersection, WriteMode::update, "out_polys");
				init_cpd(OperationDiveStrategy::zoom, 1000.0, 200.0);

			}

			ViewPair op(ViewPair& in_view_pair) override {

				ViewPair out_view;

				out_view.vector_views["out_polys"] = boost::get<GeoVector<Boost_Polygon_2>>(in_view_pair.vector_views["polys"]);
				GeoVector<Boost_Polygon_2>& gvec = boost::get<GeoVector<Boost_Polygon_2>>(out_view.vector_views["out_polys"]);

				RasterPixelsStitcher rps = RasterPixelsStitcher(in_view_pair.raster_views["proximity"]);

				for (auto& gwa : gvec.geometries_container) {
					double stitched_sum = rps.readPolygonPixels(gwa.get_definition(), RasterPixelsStitcherStartegy::contours);
					double confidence = stitched_sum / bg::area(gwa.get_definition().outer());
					gwa.set_double_attribute("confidence", confidence);
				}
				return out_view;

			}

		};

	}
}