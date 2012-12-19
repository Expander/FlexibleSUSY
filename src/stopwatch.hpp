
#include <ctime>

class Stopwatch {
public:
   Stopwatch();
   ~Stopwatch();

   void start();
   void stop();
   double get_clicks();
   double get_time_in_seconds();

private:
   clock_t time;
};
