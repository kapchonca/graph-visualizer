#include "h/graphMath.h"

double GraphMath::EuclideanDistance(Vertex* v, Vertex* u, float power) {
  return std::pow((std::pow(v->x - u->x, 2)) + pow(v->y - u->y, 2), power);
}

double GraphMath::OneVariableDerivative(Vertex* parameter, char with_respect) {
  double coordinate_p = (with_respect == 'x') ? parameter->x : parameter->y;
  double derivative = 0.0;
  for (auto neighbor : parameter->neighboorhood) {
    double coordinate_n = (with_respect == 'x') ? neighbor->x : neighbor->y;
    derivative += ((coordinate_p - coordinate_n) -
                   (kEdgeLen * parameter->distances[neighbor]) *
                       (coordinate_p - coordinate_n) /
                       EuclideanDistance(parameter, neighbor, 0.5)) /
                  parameter->distances[neighbor];
  }
  return derivative;
}

double GraphMath::CalculateDelta(Vertex* parameter) {
  double first_derivative = OneVariableDerivative(parameter, 'x');
  double second_derivative = OneVariableDerivative(parameter, 'y');
  double delta =
      std::sqrt(std::pow(first_derivative, 2) + std::pow(second_derivative, 2));
  return delta;
}

double GraphMath::TwoVariablesDerivative(Vertex* parameter, char with_respect1,
                                         char with_respect2) {
  double derivative = 0.0;
  double coordinate_p1 = (with_respect1 == 'x') ? parameter->x : parameter->y;
  double coordinate_p2 = (with_respect2 == 'x') ? parameter->x : parameter->y;
  for (auto n : parameter->neighboorhood) {
    double coordinate_n1 = (with_respect1 == 'x') ? n->x : n->y;
    double coordinate_n2 = (with_respect2 == 'x') ? n->x : n->y;
    double intermed_value =
        (kEdgeLen * parameter->distances[n] * (coordinate_p1 - coordinate_n1) *
         (coordinate_p2 - coordinate_n2)) /
        EuclideanDistance(parameter, n, 1.5);
    if (!(with_respect1 == 'x' && with_respect2 == 'y')) {
      intermed_value = 1 - intermed_value;
    }
    derivative += intermed_value / parameter->distances[n];
  }
  return derivative;
}

void GraphMath::SolveLinearEquations(Vertex* p) {

  // Coefficients of the linear equations
  double a_1 = TwoVariablesDerivative(p, 'y', 'y');
  double b_1 = TwoVariablesDerivative(p, 'x', 'y');
  double c_1 = -OneVariableDerivative(p, 'x');

  double a_2 = TwoVariablesDerivative(p, 'x', 'y');
  double b_2 = TwoVariablesDerivative(p, 'x', 'x');
  double c_2 = -OneVariableDerivative(p, 'y');

  double delta_x = 0;
  double delta_y = 0;

  if (a_1 != 0) {  // Normalize first coeff
    b_1 /= a_1;
    c_1 /= a_1;
    a_1 = 1.0;

    if (a_2 != 0) {  // Normalize first coeff
      b_2 /= a_2;
      c_2 /= a_2;
      a_2 = 1.0;

      // Subtract the 2nd equation from the 1st equation
      a_1 = 0.0;
      b_1 -= b_2;
      c_1 -= c_2;

      if (b_1 != 0) {
        delta_y = c_1 / b_1;
        delta_x = c_2 - delta_y * b_2;
      }

    } else if (b_2 != 0) {
      delta_y = c_2 / b_2;
      delta_x = (c_1 - delta_y * b_1) / a_1;
    }

  } else if (b_1 != 0) {
    delta_y = c_1 / b_1;

    if (a_2 != 0) {
      delta_x = (c_2 - delta_y * b_2) / a_2;
    }
  }

  p->x += delta_x;
  p->y += delta_y;
}
