#pragma once
#ifdef LX_SPATIAL_OPS_STATIC
#define LX_SPATIAL_OPS_API  
#else

#ifdef LxSpatialOps_EXPORTS
#define LX_SPATIAL_OPS_API __declspec(dllexport)
#else
#define LX_SPATIAL_OPS_API __declspec(dllimport)
#endif

#endif