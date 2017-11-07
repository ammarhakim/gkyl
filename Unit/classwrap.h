#ifndef CLASS_WRAP
#define  CLASS_WRAP

class Particle {
  public:
    Particle(double x) : x(x) {}
    double getx() { return x; }
  private:
    double x;
};

#endif // CLASS_WRAP
