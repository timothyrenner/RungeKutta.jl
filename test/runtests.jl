using RungeKutta
using Base.Test

# TEST 1 - Correct answer.
x_truth_test1 = [1.0 0.5416666666666; 0.0 -0.8333333333333333];
t_truth_test1 = [0.0 1.0];

f_test1 = [(t,x) -> x[2], (t,x) -> -x[1]];
x0_test1 = [1.0, 0.0];
t0_test1 = 0.0;
h_test1 = 1.0;
n_test1 = 1;

t_answer_test1, x_answer_test1 =
    rk4f(f_test1, x0_test1, t0_test1, h_test1, n_test1);

@test_approx_eq t_answer_test1 t_truth_test1
@test_approx_eq x_answer_test1 x_truth_test1

# TEST 2 - Invalid argument dimension.
f_test2 = [(t,x) -> x[2], (t,x) -> -x[1]];
x0_test2 = [1.0, 0.0, -1.0];
t0_test2 = 0.0;
h_test2 = 1.0;
n_test2 = 1;

@test_throws ArgumentError rk4f(
    f_test2, x0_test2, t0_test2, h_test2, n_test2);

# TEST 3 - Invalid number of iterations.
f_test3 = f_test2;
x0_test3 = x0_test2;
t0_test3 = t0_test2;
h_test3 = h_test2;
n_test3 = -4;

@test_throws ArgumentError rk4f(
    f_test3, x0_test3, t0_test3, h_test3, n_test3);

# TEST 4 - Invalid step size.
f_test4 = f_test3;
x0_test4 = x0_test3;
t0_test4 = t0_test3;
h_test4 = -0.5;
n_test4 = n_test2;

@test_throws ArgumentError rk4f(
    f_test4, x0_test4, t0_test4, h_test4, n_test4);
