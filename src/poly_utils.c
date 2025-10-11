#include "poly_utils.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Poly create_poly(void) {
  Poly p;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    p.coeffs[i] = 0.0;
  }
  p.degree = 0;
  return p;
}

double positive_fmod(double x, double m) {
  assert(m > 0.0);
  double r = fmod(x, m);
  if (r < 0.0)
    r += m;
  return r;
}

int64_t poly_degree(Poly p) {
  for (int64_t i = MAX_POLY_DEGREE - 1; i >= 0; i--) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      return i;
    }
  }
  return 0;
}

double get_coeff(Poly p, int64_t degree) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return 0.0;
  }
  return p.coeffs[degree];
}

void set_coeff(Poly *p, int64_t degree, double value) {
  if (degree >= MAX_POLY_DEGREE || degree < 0) {
    return;
  }
  p->coeffs[degree] = value;
  if (degree == p->degree && fabs(value) <= 1e-9) {
    while (p->degree > 0 && fabs(p->coeffs[p->degree]) <= 1e-9) {
      p->degree--;
    }
  } else if (degree > p->degree && fabs(value) > 1e-9) {
    p->degree = degree;
  }
}

Poly coeff_mod(Poly p, double modulus) {
  Poly out = create_poly();
  for (int i = 0; i <= p.degree; i++) {
    if (fabs(p.coeffs[i]) > 1e-9) {
      double rounded = round(p.coeffs[i]);
      double m = positive_fmod(rounded, modulus);
      out.coeffs[i] = m;
    }
  }
  out.degree = p.degree;
  return out;
}

Poly poly_add(Poly a, Poly b) {
  Poly sum = create_poly();
  int max_degree = (a.degree > b.degree) ? a.degree : b.degree;
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    sum.coeffs[i] = a.coeffs[i] + b.coeffs[i];
  }
  int64_t deg = max_degree;
  while (deg > 0 && fabs(sum.coeffs[deg]) < 1e-9) deg--;
  sum.degree = deg;
  return sum;
}

Poly poly_mul_scalar(Poly p, double scalar) {
  Poly res = create_poly();
  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    res.coeffs[i] = p.coeffs[i] * scalar;
    if (fabs(res.coeffs[i]) > 1e-9) {
        res.degree = i;
    }
  }
  return res;
}

Poly poly_mul(Poly a, Poly b) {
  Poly res = create_poly();

  int nonzero_deg[MAX_POLY_DEGREE];
  size_t nz_count = 0;
  for (int i = 0; i <= b.degree; i++) {
      if (fabs(b.coeffs[i]) > 1e-9)
          nonzero_deg[nz_count++] = i;
  }
  int max_res_degree = 0;
  for (int i = 0; i <= a.degree; i++) {
    if (fabs(a.coeffs[i]) > 1e-9) {
      for (size_t j = 0; j < nz_count; j++) {
        int ind = nonzero_deg[j];
        assert(i + ind < MAX_POLY_DEGREE);
        res.coeffs[i + ind] += a.coeffs[i] * b.coeffs[ind];
        if (i + ind > max_res_degree && fabs(res.coeffs[i + ind]) > 1e-9)
            max_res_degree = i + ind;
      }
    }
  }
  res.degree = max_res_degree;
  return res;
}

void poly_divmod(Poly num, Poly den, Poly *quot, Poly *rem) {
  // In our case `den` should always be (x^n + 1)
  assert(poly_degree(den) > 0 || fabs(get_coeff(den, 0)) > 1e-9);

  size_t ndeg = poly_degree(num);
  size_t ddeg = poly_degree(den);

  *quot = create_poly();
  *rem = num;

  if (ndeg < ddeg) {
    return;
  }

  int nonzero_deg[MAX_POLY_DEGREE];
  size_t nz_count = 0;
  for (int i = 0; i <= ddeg; i++) {
      if (fabs(den.coeffs[i]) > 1e-9)
          nonzero_deg[nz_count++] = i;
  }

  double d_lead = get_coeff(den, ddeg);
  assert(fabs(d_lead) > 1e-9);
  int max_rem_degree = rem->degree;
  int max_quot_degree = quot->degree;
  for (int64_t k = ndeg - ddeg; k >= 0; --k) {
    int64_t target_deg = ddeg + k;
    double r_coeff = get_coeff(*rem, target_deg);
    double coeff = trunc(round(r_coeff) / round(d_lead));
    quot->coeffs[k] += coeff;
    if (k > max_quot_degree && fabs(quot->coeffs[k]) > 1e-9)
        max_quot_degree = k;

    for (size_t j = 0; j < nz_count; j++) {
        int i = nonzero_deg[j];
        rem->coeffs[i + k] -= coeff * den.coeffs[i];
        if (i + k > max_rem_degree && fabs(rem->coeffs[i + k]) > 1e-9)
            max_rem_degree = i + k;
    }
  }
  quot->degree = max_quot_degree;
  rem->degree = max_rem_degree;

  assert(poly_degree(*rem) < poly_degree(den));
}

Poly poly_round_div_scalar(Poly x, double divisor) {
  Poly out = create_poly();
  assert(fabs(divisor) > 1e-9);

  for (int i = 0; i < MAX_POLY_DEGREE; i++) {
    double v = x.coeffs[i];
    out.coeffs[i] = round(v / divisor);
  }
  int64_t deg = x.degree;
  while (deg > 0 && fabs(out.coeffs[deg]) < 1e-9) deg--;
  out.degree = deg;
  return out;
}
