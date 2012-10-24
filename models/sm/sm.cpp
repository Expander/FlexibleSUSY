
#include "sm.hpp"

#include <cassert>

StandardModel::StandardModel()
   : yu(3, 3), yd(3, 3), ye(3, 3), g(3)
{
   setPars(numStandardModelPars);
   setMu(0.0);
   setLoops(1);
   setThresholds(0);
}

StandardModel::StandardModel(const StandardModel& s)
   : yu(s.yu), yd(s.yd), ye(s.ye), g(s.g)
{
   setPars(numStandardModelPars);
   setMu(s.displayMu());
   setLoops(s.displayLoops());
   setThresholds(s.displayThresholds());
}

StandardModel::StandardModel(const DoubleMatrix& SMu, const DoubleMatrix& SMd,
                             const DoubleMatrix& SMe, const DoubleVector& g_)
   : yu(SMu), yd(SMd), ye(SMe), g(g_)
{
   setPars(numStandardModelPars);
   setMu(0.0);
   setLoops(1);
   setThresholds(0);
}

StandardModel::~StandardModel()
{
}

const StandardModel& StandardModel::operator=(const StandardModel& s)
{
   if (this == &s) return *this;
   yu = s.yu;
   yd = s.yd;
   ye = s.ye;
   g = s.g;
   setMu(s.displayMu());
   setLoops(s.displayLoops());
   setThresholds(s.displayThresholds());
   return *this;
}

void StandardModel::setGaugeCoupling(int i, double f)
{
   g(i) = f;
}

void StandardModel::setAllGauge(const DoubleVector& v)
{
   assert(v.displayStart() == 1 && v.displayEnd() == 3);
   g = v;
}

DoubleVector StandardModel::displayGauge() const
{
   return g;
}

double StandardModel::displayGaugeCoupling(int i) const
{
   return g.display(i);
}

void StandardModel::setYukawaElement(yukawa k, int i, int j, double f)
{
   switch (k) {
   case YU:
      yu(i, j) = f;
      break;
   case YD:
      yd(i, j) = f;
      break;
   case YE:
      ye(i, j) = f;
      break;
   default:
      assert(false && "StandardModel::setYukawaElement called with illegal k");
      break;
   }
}

void StandardModel::setYukawaMatrix(yukawa k, const DoubleMatrix& m)
{
   switch (k) {
   case YU:
      yu = m;
      break;
   case YD:
      yd = m;
      break;
   case YE:
      ye = m;
      break;
   default:
      assert(false && "StandardModel::setYukawaMatrix called with illegal k");
      break;
   }
}

double StandardModel::displayYukawaElement(yukawa k, int i, int j) const
{
   switch (k) {
   case YU:
      return yu.display(i, j);
      break;
   case YD:
      return yd.display(i, j);
      break;
   case YE:
      return ye.display(i, j);
      break;
   default:
      assert(false && "StandardModel::displayYukawaElement called with illegal k");
      break;
   }
   return 0.0;
}

DoubleMatrix StandardModel::displayYukawaMatrix(yukawa k) const
{
   switch (k) {
   case YU:
      return yu;
      break;
   case YD:
      return yd;
      break;
   case YE:
      return ye;
      break;
   default:
      assert(false && "StandardModel::displayYukawaMatrix called with illegal k");
      break;
   }
}

//Peter:: edited to include new Essm parameters

const DoubleVector StandardModel::display() const
{
   DoubleVector y(numStandardModelPars);
   int i, j, k = 0;
   for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++) {
         k++;
         y(k) = yu.display(i, j);
         y(k + 9) = yd.display(i, j);
         y(k + 18) = ye.display(i, j);
      }
   k = 27;
   for (i = 1; i <= 3; i++) {
      k++;
      y(k) = g.display(i);
   }

   return y;
}
//Peter:: edited to include new Essm parameters

void StandardModel::set(const DoubleVector& y)
{
   int i, j, k = 0;
   for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++) {
         k++;
         yu(i, j) = y.display(k);
         yd(i, j) = y.display(k + 9);
         ye(i, j) = y.display(k + 18);
      }
   k = 27;
   for (i = 1; i <= 3; i++) {
      k++;
      g(i) = y.display(k);
   }
}

std::ostream& operator <<(std::ostream& left, const StandardModel& s)
{
   left << "SM parameters at Q: " << s.displayMu()
        << '\n'
        << " Y^U" << s.displayYukawaMatrix(StandardModel::YU)
        << " Y^D" << s.displayYukawaMatrix(StandardModel::YD)
        << " Y^E" << s.displayYukawaMatrix(StandardModel::YE)
        << '\n'
        << " g1: " << s.displayGaugeCoupling(1)
        << " g2: " << s.displayGaugeCoupling(2)
        << " g3: " << s.displayGaugeCoupling(3)
        << '\n'
        << " thresholds: " << s.displayThresholds()
        << " #loops: " << s.displayLoops() << '\n';
   return left;
}

// Outputs derivatives (DRbar scheme) in the form of ds
StandardModel StandardModel::calcBeta() const
{
   static const double oneO16Pisq = 1.0 / (16.0 * PI * PI);
   DoubleMatrix dyu(3, 3), dyd(3, 3), dye(3, 3);
   DoubleVector dg(3);

   dyu(3, 3) = oneO16Pisq * yu.display(3, 3) * (
      -17.0 / 20.0 * sqr(displayGaugeCoupling(1))
      - 9.0 / 4.0 * sqr(displayGaugeCoupling(2))
      - 8.0 * sqr(displayGaugeCoupling(3))
      + 4.5 * sqr(yu.display(3, 3))
      + 1.5 * sqr(yd.display(3, 3))
      + sqr(ye.display(3, 3)));

   dyd(3, 3) = oneO16Pisq * yd.display(3, 3) * (
      -0.25 * sqr(displayGaugeCoupling(1))
      - 9.0 / 4.0 * sqr(displayGaugeCoupling(2))
      - 8.0 * sqr(displayGaugeCoupling(3))
      + 1.5 * sqr(yu.display(3, 3))
      + 4.5 * sqr(yd.display(3, 3))
      + sqr(ye.display(3, 3)));

   dye(3, 3) = oneO16Pisq * ye.display(3, 3) * (
      -9.0 / 4.0 * sqr(displayGaugeCoupling(1))
      -9.0 / 4.0 * sqr(displayGaugeCoupling(2))
      + 3.0 * sqr(yu.display(3, 3))
      + 3.0 * sqr(yd.display(3, 3))
      + 2.5 * sqr(ye.display(3, 3)));

   dg(1) = oneO16Pisq * std::pow(displayGaugeCoupling(1), 3) * (41.0 / 10.0);
   dg(2) = oneO16Pisq * std::pow(displayGaugeCoupling(2), 3) * (-19.0 / 6.0);
   dg(3) = oneO16Pisq * std::pow(displayGaugeCoupling(3), 3) * (-7.0);

   return StandardModel(dyu, dyd, dye, dg);
}

DoubleVector StandardModel::beta() const
{
   return calcBeta().display();
}
