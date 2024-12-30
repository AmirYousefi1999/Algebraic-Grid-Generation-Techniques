# Algebraic-Grid-Generation-Techniques
MATLAB Code for Algebraic Grid Generation with Partial Derivatives and Jacobian Calculation Below is a MATLAB code that generates an algebraic grid for a 
rectangular domain, calculates partial derivatives, computes the Jacobian at each grid point, and evaluates the error compared to exact analytical values.


       clc;
      clear;
      H1 = 2; H2 = 4; L = 4;
      Im = 16; Jm = 12;
      dkes = 1 / (Im - 1); deta = 1 / (Jm - 1);

     for i = 1:Im
        for j = 1:Jm
      kes = dkes * (i - 1); eta = deta * (j - 1);
      x(i, j) = kes * L;
      y(i, j) = (H1 + (H2 - H1) * kes) * eta;
        end
    end

    for i = 1:Im
      for j = 1:Jm
        kx_exact(i, j) = 1 / L;
        ky_exact(i, j) = 0.0;
        ex_exact(i, j) = -(H2 - H1) / L * y(i, j) / (H1 + (H2 - H1) / L * x(i, j))^2;
        ey_exact(i, j) = 1 / (H1 + (H2 - H1) / L * x(i, j));

          if i == 1
            xk = (-3 * x(i, j) + 4 * x(i + 1, j) - x(i + 2, j)) / (2 * dkes);
            yk = (-3 * y(i, j) + 4 * y(i + 1, j) - y(i + 2, j)) / (2 * dkes);
        elseif i == Im
            xk = (x(i - 2, j) - 4 * x(i - 1, j) + 3 * x(i, j)) / (2 * dkes);
            yk = (y(i - 2, j) - 4 * y(i - 1, j) + 3 * y(i, j)) / (2 * dkes);
        else
            xk = (x(i + 1, j) - x(i - 1, j)) / (2 * dkes);
            yk = (y(i + 1, j) - y(i - 1, j)) / (2 * dkes);
        end
        if j == 1
            xe = (-3 * x(i, j) + 4 * x(i, j + 1) - x(i, j + 2)) / (2 * deta);
            ye = (-3 * y(i, j) + 4 * y(i, j + 1) - y(i, j + 2)) / (2 * deta);
        elseif j == Jm
            xe = (x(i, j - 2) - 4 * x(i, j - 1) + 3 * x(i, j)) / (2 * deta);
            ye = (y(i, j - 2) - 4 * y(i, j - 1) + 3 * y(i, j)) / (2 * deta);
        else
            xe = (x(i, j + 1) - x(i, j - 1)) / (2 * deta);
            ye = (y(i, j + 1) - y(i, j - 1)) / (2 * deta);
        end
        Jac = 1 / (xk * ye - yk * xe);
        kx(i, j) = Jac * ye;
        ky(i, j) = -Jac * xe;
        ex(i, j) = -Jac * yk;
        ey(i, j) = Jac * xk;
        error_kx(i, j) = kx(i, j) - kx_exact(i, j);
        error_ky(i, j) = ky(i, j) - ky_exact(i, j);
        error_ex(i, j) = ex(i, j) - ex_exact(i, j);
        error_ey(i, j) = ey(i, j) - ey_exact(i, j);

    end
end

surface(x, y, y * 0);

