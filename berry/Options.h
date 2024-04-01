#pragma once

/* Enable bounds checking */
#define BRY_ENABLE_BOUNDS_CHECK

/* Enable inline */
#define BRY_ENABLE_INL

/* Enable asserts */
#define BRY_ASSERTS

/* Enable debug tools */
#define BRY_DEBUG_TOOLS

/* Enable logging in color */
#define BRY_LOG_COLOR


#define EIGEN_FFTW_DEFAULT



#ifdef BRY_ENABLE_INL
    #define BRY_INL inline
#elif
    #define BRY_INL
#endif