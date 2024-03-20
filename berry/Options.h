#pragma once

/* Enable bounds checking */
#define BRY_ENABLE_BOUNDS_CHECK

/* Enable inline */
#define BRY_ENABLE_INL







#ifdef BRY_ENABLE_INL
    #define BRY_INL inline
#elif
    #define BRY_INL
#endif