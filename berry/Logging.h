#pragma once

#include <iostream>

#include "Options.h"

#ifdef BRY_LOG_COLOR
    #define BRY_LOG_WHITE(msg) "\033[0;37m" << msg << "\033[0m"
    #define BRY_LOG_BWHITE(msg) "\033[1;37m" << msg << "\033[0m"
    #define BRY_LOG_GREEN(msg) "\033[0;32m" << msg << "\033[0m"
    #define BRY_LOG_BGREEN(msg) "\033[1;32m" << msg << "\033[0m"
    #define BRY_LOG_YELLOW(msg) "\033[0;33m" << msg << "\033[0m"
    #define BRY_LOG_BYELLOW(msg) "\033[1;33m" << msg << "\033[0m"
    #define BRY_LOG_RED(msg) "\033[0;31m" << msg << "\033[0m"
    #define BRY_LOG_BRED(msg) "\033[1;31m" << msg << "\033[0m"
    #define BRY_LOG_CYAN(msg) "\033[0;36m" << msg << "\033[0m"
    #define BRY_LOG_BCYAN(msg) "\033[1;36m" << msg << "\033[0m"
    #define BRY_LOG_CYAN(msg) "\033[0;36m" << msg << "\033[0m"
    #define BRY_LOG_BCYAN(msg) "\033[1;36m" << msg << "\033[0m"
#else
    #define BRY_LOG_WHITE(msg) msg
    #define BRY_LOG_BWHITE(msg) msg
    #define BRY_LOG_GREEN(msg) msg
    #define BRY_LOG_BGREEN(msg) msg
    #define BRY_LOG_YELLOW(msg) msg
    #define BRY_LOG_BYELLOW(msg) msg
    #define BRY_LOG_RED(msg) msg
    #define BRY_LOG_BRED(msg) msg
    #define BRY_LOG_CYAN(msg) msg
    #define BRY_LOG_BCYAN(msg) msg
    #define BRY_LOG_CYAN(msg) msg
    #define BRY_LOG_BCYAN(msg) msg
#endif

#ifdef BRY_EXCLUDE_LOGS
    #define LOG(msg) 
    #define PRINT(msg) 
    #define PRINT_VEC2(msg, vec2)
    #define PRINT_VEC3(msg, vec3)
    #define PRINT_NAMED(name, msg) 
    #define DEBUG(msg)
    #define INFO(msg) 
    // Keep error message
    #define ERROR(msg) std::cerr << "\033[1;31m >[ERR] ERROR ("<< __func__ << "): \033[1;31m" << msg << "\033[0m \n"
    #define WARN(msg) 
    #define NEW_LINE 
#else
    #define LOG(msg) std::cout << "\033[1;36m >[LOG]\033[0;37m " << msg << "\033[0m \n"
    #define PRINT(msg) std::cout << "\033[0;37m" << msg << "\033[0m \n"
    #define PRINT_VEC2(msg, vec2) std::cout << "\033[0;37m" << msg << " (" << vec2[0] << ", " << vec2[1] << ")\033[0m \n"
    #define PRINT_VEC3(msg, vec3) std::cout << "\033[0;37m" << msg << " (" << vec3[0] << ", " << vec3[1] << ", " << vec3[2] << ")\033[0m \n"
    #define PRINT_NAMED(name, msg) std::cout << "        "<< BRY_LOG_BGREEN(name) << ": " << BRY_LOG_WHITE(msg)
    #define DEBUG(msg) std::cout << BRY_LOG_BCYAN(" >[DBG] ("<< __func__ << "): ") << BRY_LOG_WHITE(msg) << "\n"
    #define INFO(msg) std::cout << BRY_LOG_BWHITE(" >[IFO]: ") << BRY_LOG_WHITE(msg) << " \n"
    #define ERROR(msg) std::cerr << BRY_LOG_BRED(" >[ERR] ERROR ("<< __func__ << "): ") << BRY_LOG_WHITE(msg) << "\n"
    #define WARN(msg) std::cout << BRY_LOG_BYELLOW(" >[WRN] WARNING ("<< __func__ << "): ") << BRY_LOG_WHITE(msg) << "\n"
    #define NEW_LINE std::cout<<"\n"
#endif

#ifdef BRY_ASSERTS
    #define ASSERT(condition, msg) {if (!(condition)) {ERROR("[Assert fail] " << msg); exit(1);}}
#endif

#ifdef BRY_DEBUG_TOOLS
    #define PAUSE std::cout << "<paused>" << "\n"; \
            std::cin.get()
#endif

