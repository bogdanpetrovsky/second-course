namespace second_course;

public static class FunctionHelper
{
    // double[] X1(double t) { return new [] { Math.Cos(t), Math.Sin(t) }; }
    // double[] X2(double t) { return new [] { 2 * Math.Cos(t), 2 * Math.Sin(t) }; }
    // double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), Math.Cos(t) }; }
    // double[] Der1X2(double t) { return new [] { -2 * Math.Sin(t), 2 * Math.Cos(t) }; }
    // double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -1 * Math.Sin(t) }; }
    // double[] Der2X2(double t) { return new [] { -2 * Math.Cos(t), -2 * Math.Sin(t) }; }
    // double[] Der3X1(double t) { return new [] { Math.Sin(t), -1 * Math.Cos(t) }; }
    // double[] Der3X2(double t) { return new [] { 2 * Math.Sin(t), -2 * Math.Cos(t) }; }
    public static double[] X1(double t) { return new [] { Math.Cos(t), 0.5 * Math.Sin(t) }; }
    public static double[] X2(double t) { return new [] { 1.5 * Math.Cos(t), Math.Sin(t) }; }
    public static double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), 0.5 * Math.Cos(t) }; }
    public static double[] Der1X2(double t) { return new [] { -1.5 * Math.Sin(t), Math.Cos(t) }; }
    public static double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -0.5 * Math.Sin(t) }; }
    public static double[] Der2X2(double t) { return new [] { -1.5 * Math.Cos(t), -1 * Math.Sin(t) }; }
    public static double[] Der3X1(double t) { return new [] { Math.Sin(t), -0.5 * Math.Cos(t) }; }
    public static double[] Der3X2(double t) { return new [] { 1.5 * Math.Sin(t), -1 * Math.Cos(t) }; }
    public static double Eps = 0.000001;
    public static double ErrorEps = 0.01;
    
    public static double[] VGamma2(double t)
    {
        return new []
        {
            Der1X2(t)[1] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0),
            -1 * Der1X2(t)[0] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0)
        };
    }
    public static double[] VGamma1(double t)
    {
        return new []
        {
            Der1X1(t)[1] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0),
            -1 * Der1X1(t)[0] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0)
        };
    }
    
    public static double GetEuclideanDistance(double x21, double x11, double x22, double x12)
    {
        return Math.Sqrt(Math.Pow(x21 - x11, 2) + Math.Pow(x22 - x12, 2));
    }
    
    public static double[] Gauss (double[,] a, int N) {
        int n = N;
        int m = N;
 
        int[] where = new int[m];
        for (int i=0; i<m; i++) { where[i] = -1; }
    
        for (int col=0, row=0; col<m && row<n; ++col) {
            int sel = row;
            for (int i=row; i<n; ++i)
                if (Math.Abs(a[i, col]) > Math.Abs(a[sel, col]))
                    sel = i;
            if (Math.Abs (a[sel, col]) < Eps)
                continue;
            for (int i = col; i <= m; ++i)
                (a[sel, i], a[row, i]) = (a[row, i], a[sel, i]);
            where[col] = row;
 
            for (int i=0; i<n; ++i)
                if (i != row) {
                    double c = a[i, col] / a[row, col];
                    for (int j=col; j<=m; ++j)
                        a[i, j] -= a[row, j] * c;
                }
            ++row;
        }
    
        // for (int i = 0; i < m; i++) { Console.Write(where[i] + " "); } Console.WriteLine();
    
        double[] ans = new double[m];
        for (int i=0; i<m; i++) { ans[i] = 0; }
    
        for (int i=0; i<m; ++i)
            if (where[i] != -1)
                ans[i] = a[where[i], m] / a[where[i], i];
        for (int i=0; i<n; ++i) {
            double sum = 0;
            for (int j=0; j<m; ++j)
                sum += ans[j] * a[i, j];
            if (Math.Abs(sum - a[i, m]) > Eps)
                return ans;
        }
 
        for (int i=0; i<m; ++i)
            if (where[i] == -1)
                return ans;
        return ans;
    }
    
    public static double Rj(double t, double tj, double n)
    {
        double result = 0;
        
        for (int m = 1; m < n; m++)
        {
            result += 1.0 / m * Math.Cos(m * (t - tj)) + 1 / (2 * n) * Math.Cos(n * (t - tj));

        }

        return - result / n;
    }

    public static double GetApproximatedU(double x1, double x2, double[] uValues, int N, bool isReversedProblem)
    {
        double s1 = 0, s2 = 0;
        for (int i = 0; i < 2*N; i++)
        {
            double t = i * Math.PI / N;
            if (isReversedProblem)
            {
                s1 = s1 + uValues[i] * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0) *
                    (((x1 - X1(t)[0]) * VGamma1(t)[0] + (x2 - X1(t)[1]) * VGamma1(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1]), 2));
                s2 = s2 + uValues[i+2*N] * Math.Log(1 / GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1])) * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0);
            }
            else
            {
                s1 = s1 + uValues[i] * Math.Log(1 / GetEuclideanDistance(x1, X1(t)[0], x2, X1(t)[1])) * GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0);
                s2 = s2 + uValues[i+2*N] * GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0) *
                      (((x1 - X2(t)[0]) * VGamma2(t)[0] + (x2 - X2(t)[1]) * VGamma2(t)[1]) / Math.Pow(GetEuclideanDistance(x1, X2(t)[0], x2, X2(t)[1]), 2));
            }
        }

        return s1 / (2 * N) + s2 / (2 * N);
    }

    public static double[] GetApproximatedUOnGamma1(double[] uValues, int N)
    {
        double[] result = new double[2*N];
        
        for (int k = 0; k < 2*N; k++)
        {
            double tk = k * Math.PI / N;
            double s1 = 0, s2 = 0;
            for (int i = 0; i < 2*N; i++)
            {
                double ti = i * Math.PI / N;
                s2 = s2 + uValues[i+2*N] * GetEuclideanDistance(Der1X2(ti)[0], 0, Der1X2(ti)[1], 0) *
                    ((X1(tk)[0] - X2(ti)[0]) * VGamma2(ti)[0] + (X1(tk)[1] - X2(ti)[1]) * VGamma2(ti)[1]) / Math.Pow(GetEuclideanDistance(X1(tk)[0], X2(ti)[0], X1(tk)[1], X2(ti)[1]), 2);

                double l = 0;
                if (Math.Abs(ti - tk) < Eps)
                {
                    l = Math.Log(Math.Pow(GetEuclideanDistance(Der1X1(ti)[0], 0, Der1X1(ti)[1], 0), 2));
                }
                else
                {
                    l = Math.Log(Math.Pow(GetEuclideanDistance(X1(tk)[0], X1(ti)[0], X1(tk)[1], X1(ti)[1]), 2) / (4 * Math.Pow(Math.Sin((tk - ti)/2), 2)));
                }

                s1 = s1 + uValues[i] * GetEuclideanDistance(Der1X1(ti)[0], 0, Der1X1(ti)[1], 0) * (-l/(4*N) - Rj(ti,tk, N)/2);
            }

            result[k] = s1 + s2 / (2 * N);
        }


        return result;
    }

    public static double HWaved(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < Eps)
        {
            return -1.0 / 6
                   +
                   (Der1X1(t1)[0] * Der3X1(t1)[0] + Der1X1(t2)[1] * Der3X1(t1)[1]) /
                   (3 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2))
                   +
                   Math.Pow(GetEuclideanDistance(Der2X1(t1)[0], 0, Der2X1(t1)[1], 0), 2) /
                   (2 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2))
                   -
                   Math.Pow(Der1X1(t1)[0] * Der2X1(t1)[0] + Der1X1(t2)[1] * Der2X1(t1)[1], 2) /
                   Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 4);
        }
        else
        {
            double denominator1 = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 4);
            double numerator1 = (Der1X1(t1)[0] * (X1(t2)[0] - X1(t1)[0]) + Der1X1(t1)[1] * (X1(t2)[1] - X1(t1)[1])) *
                                (Der1X1(t2)[0] * (X1(t1)[0] - X1(t2)[0]) + Der1X1(t2)[1] * (X1(t1)[1] - X1(t2)[1]));
            double denominator2 = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 2);
            double numerator2 = (Der1X1(t1)[0] * Der1X1(t2)[0] + Der1X1(t1)[1] * Der1X1(t2)[1]);
            
            double denominator3 = Math.Pow(Math.Sin((t1-t2)/2), 2);
            return 4 * numerator1 / denominator1 - 2 * numerator2 / denominator2 - 1 / (2 * denominator3);
        }
    }

    public static double T1J(double t, double tj, int n)
    {
        double sum = 0;
        for (int m = 1; m < n; m++)
        {
            sum += m * Math.Cos(m * (t - tj));
        }
        
        return -sum / n - Math.Cos(n * (t - tj)) / 2;
    }

    public static double T2J(double t1, double t2)
    {
        double kernelNumerator = (X2(t2)[0] - X1(t1)[0]) * VGamma1(t1)[0] + (X2(t2)[1] - X1(t1)[1]) * VGamma1(t1)[1];
        double kernelDenominator = Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 2);
        
        return kernelNumerator * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0) /
               kernelDenominator;
    }

    public static double[] GetApproximatedDerUOnGamma1(double[] uValues, int N)
    {
        double[] result = new double[2*N];
        
        for (int k = 0; k < 2*N; k++)
        {
            double tk = k * Math.PI / N;
            double s1 = 0, s2 = 0;
            for (int i = 0; i < 2*N; i++)
            {
                double ti = i * Math.PI / N;
                s2 = s2 + uValues[i + 2 * N] * T2J(tk, ti);
                s1 = s1 + (uValues[i] * T1J(tk, ti, N) + uValues[i] * HWaved(tk, ti)/(2*N)) /
                    GetEuclideanDistance(Der1X1(ti)[0], 0, Der1X1(ti)[1], 0);;
            }

            result[k] = s1 + s2 / (2 * N);
        }


        return result;
    }

    public static double F1(double x)
    {
        return Math.Pow(X1(x)[0], 2) - Math.Pow(X1(x)[1], 2);
    }

    public static double G2(double x)
    {
        return 2 * X2(x)[0] * VGamma2(x)[0] - 2 * X2(x)[1] * VGamma2(x)[1];
    }

    public static double du_G1(double x)
    {
        return 2 * X1(x)[0] * VGamma1(x)[0] - 2 * X1(x)[1] * VGamma1(x)[1] ;

    }
    
    public static double GetMaxError(double[] f1Value, double[] f2Value)
    {
        double maxError = 0;

        for (int i = 0; i < f1Value.Length; i++)
        {
            maxError = Math.Max(Math.Abs(f1Value[i] - f2Value[i]), maxError);
        }

        return maxError;
    }
}