
#include <tma.h>

using namespace tma;

int main(int argc, char**argv)
{
  // interval cell
  message("interval");
  {
    interval C;
    tma_assert(C.dim()  == 1);
    tma_assert(C.num(0) == 2);
    tma_assert(C.num(1) == 1);

    interval::reference R;
    tma_assert(R.x(0)[0] == 0.0);
    tma_assert(R.x(1)[0] == 1.0);

    // Interval with two subintervals
    topology<interval> T(2, 3);
    T(0)[0] = 0;
    T(0)[1] = 1;
    T(1)[0] = 1;
    T(1)[1] = 2;
    T.dump();
    geometry<1> G(3);
    G(0)[0] = 0.0;
    G(1)[0] = 0.5;
    G(2)[0] = 1.0;
    G.dump();

    // Test mesh constructor
    mesh<interval, 1> (2, 3);
  }

  // P1 on interval cell
  message("P1::interval");
  {
    P1::interval E;
    // Check element dimension
    tma_assert(E.dim()   == 2);
    // Check element node coordinates
    tma_assert(E.x(0)[0] == 0.0);
    tma_assert(E.x(1)[0] == 1.0);
    // Check basis functions
    real v[2];
    E(E.x(0), v);
    tma_assert(v[0]  == 1.0);
    tma_assert(v[1]  == 0.0);
    E(E.x(1), v);
    tma_assert(v[0]  == 0.0);
    tma_assert(v[1]  == 1.0);
  }

  // triangle cell
  message("triangle");
  {
    triangle C;
    tma_assert(C.dim()  == 2);
    tma_assert(C.num(0) == 3);
    tma_assert(C.num(1) == 3);
    tma_assert(C.num(2) == 1);


    triangle::reference R;
    tma_assert(R.x(0)[0] == 0.0);
    tma_assert(R.x(0)[1] == 0.0);
    tma_assert(R.x(1)[0] == 1.0);
    tma_assert(R.x(1)[1] == 0.0);
    tma_assert(R.x(2)[0] == 0.0);
    tma_assert(R.x(2)[1] == 1.0);

    // A square with two triangles would use this topology
    topology<triangle> T(2, 4);
    T(0)[0] = 0;
    T(0)[1] = 1;
    T(0)[2] = 2;
    T(1)[0] = 0;
    T(1)[1] = 2;
    T(1)[2] = 3;
    T.dump();

    geometry<2> G(4);
    G(0)[0] = 0.0; G(0)[1] = 0.0;
    G(1)[0] = 1.0; G(1)[1] = 0.0;
    G(2)[0] = 1.0; G(2)[1] = 1.0;
    G(3)[0] = 0.0; G(3)[1] = 1.0;
    G.dump();

    // Test mesh constructor
    mesh<triangle, 2> (2, 4);
  }

  // P1 on triangle cell
  message("P1::triangle");
  {
    P1::triangle E;
    // Check element dimension
    tma_assert(E.dim()   == 3);
    // Check element node coordinates
    tma_assert(E.x(0)[0] == 0.0);
    tma_assert(E.x(0)[1] == 0.0);
    tma_assert(E.x(1)[0] == 1.0);
    tma_assert(E.x(1)[1] == 0.0);
    tma_assert(E.x(2)[0] == 0.0);
    tma_assert(E.x(2)[1] == 1.0);
    // Check basis functions
    real v[3];
    E(E.x(0), v);
    tma_assert(v[0]  == 1.0);
    tma_assert(v[1]  == 0.0);
    tma_assert(v[2]  == 0.0);
    E(E.x(1), v);
    tma_assert(v[0]  == 0.0);
    tma_assert(v[1]  == 1.0);
    tma_assert(v[2]  == 0.0);
    E(E.x(2), v);
    tma_assert(v[0]  == 0.0);
    tma_assert(v[1]  == 0.0);
    tma_assert(v[2]  == 1.0);
  }

  // QR

  return 0;
}
