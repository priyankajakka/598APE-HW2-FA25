#include "he.h"
#include "poly_random.h"
#include "poly_utils.h"
#include "ring_utils.h"

Poly encode_plain_integer(double t, double pt) {
  Poly m = create_poly();
  double v = positive_fmod(pt, t);
  m.coeffs[0] = v;
  m.degree = 0;
  return m;
}

Ciphertext encrypt(PublicKey pk, size_t n, double q, Poly poly_mod, double t,
                   double pt) {
  Poly m = encode_plain_integer(t, pt);
  Poly scaled_m = poly_mul_scalar(m, floor(q / t));
  Poly e1 = gen_normal_poly(n, 0.0, 1.0);
  Poly e2 = gen_normal_poly(n, 0.0, 1.0);
  Poly u = gen_binary_poly(n);

  Poly bu = ring_mul_mod(pk.b, u, q, poly_mod);
  Poly bu_e1 = ring_add_mod(bu, e1, q, poly_mod);
  Poly c0 = ring_add_mod(bu_e1, scaled_m, q, poly_mod);

  Poly au = ring_mul_mod(pk.a, u, q, poly_mod);
  Poly c1 = ring_add_mod(au, e2, q, poly_mod);

  Ciphertext ct;
  ct.c0 = c0;
  ct.c1 = c1;
  return ct;
}