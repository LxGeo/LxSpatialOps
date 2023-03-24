#pragma once
#include "defs.h"
#include "base_operation.h"
#include "spatial_datasets/patchified_dataset.h"

namespace LxGeo
{
	namespace LxRasterOps
	{

		class PolygonToLineString : public BaseOperation {

		public:
			PolygonToLineString(
				const std::string& in_vec,
				const std::string& out_vec
			) {

				add_vector_input_dataset(in_vec, "in_vec");
				add_vector_output_dataset(out_vec, { "in_vec" }, ExtentsCombinationStrategy::ext_union, "out_vec");
				init_cpd(OperationDiveStrategy::same);
			};

			ViewPair op(ViewPair& in_view_pair) override {
				ViewPair out_view;
				auto& v_views = in_view_pair.vector_views;
				const GeoVector<Boost_Polygon_2>& gvec = boost::get<GeoVector<Boost_Polygon_2>>(v_views["in_vec"]);
				GeoVector<Boost_LineString_2> out_vec; 

				std::function<Geometries_with_attributes<Boost_LineString_2>(const Geometries_with_attributes<Boost_Polygon_2>&)> transformer_fn = [](const Geometries_with_attributes<Boost_Polygon_2>& in_gwa)->Geometries_with_attributes<Boost_LineString_2> {
					Geometries_with_attributes<Boost_LineString_2> out_gwa(in_gwa); //copy all attributes
					Boost_LineString_2& out_ls = out_gwa.get_definition();
					out_ls.assign(in_gwa.get_definition().outer().begin(), in_gwa.get_definition().outer().end());
					return out_gwa;
				};

				IO_DATA::transform_geovector(gvec, out_vec, transformer_fn);

				out_view.vector_views["out_vec"] = out_vec;
				return out_view;
			}
		};

	}
}