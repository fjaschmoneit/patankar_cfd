//
// Created by Peter Berg Ammundsen on 05/01/2026.
//

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
    };
};

#endif //UTIL_H
