
#include "stopwatch.hpp"

Stopwatch::Stopwatch()
   : time(0)
{
}

Stopwatch::~Stopwatch()
{
}

void Stopwatch::start()
{
   time = clock();
}

void Stopwatch::stop()
{
   time = clock() - time;
}

double Stopwatch::get_clicks()
{
   return time;
}

double Stopwatch::get_time_in_seconds()
{
   return static_cast<double>(time)/CLOCKS_PER_SEC;
}
