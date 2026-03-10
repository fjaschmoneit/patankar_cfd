#include <chrono>

#ifndef UTIL_H
#define UTIL_H

namespace util
{
    class timer {
    private:
        std::chrono::high_resolution_clock::time_point start_time;
    public:
        void start(){
            start_time = std::chrono::high_resolution_clock::now();
        }

        long long stop() {
            auto end_time = std::chrono::high_resolution_clock::now();
            return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        }

        static std::string today() {

            const std::chrono::time_point now{std::chrono::system_clock::now()};
            const std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)};

            std::string y = std::to_string(static_cast<int>(ymd.year()));
            std::string m = std::to_string(static_cast<unsigned>(ymd.month()));
            std::string d = std::to_string(static_cast<unsigned>(ymd.day()));
            return y + m + d;
        }

    };
};

#endif //UTIL_H
