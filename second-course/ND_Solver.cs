using second_course;

public class ND_Solver
{
    double K11(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < FunctionHelper.Eps)
        {
            return (FunctionHelper.Der2X1(t1)[0] * FunctionHelper.Der1X1(t1)[1] - FunctionHelper.Der2X1(t1)[1] * FunctionHelper.Der1X1(t1)[0]) /
                   (2 * Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t1)[0], 0, FunctionHelper.Der1X1(t1)[1], 0), 2));
        }

        double kernelNumerator = (FunctionHelper.X1(t2)[0] - FunctionHelper.X1(t1)[0]) * FunctionHelper.VGamma1(t1)[0] + (FunctionHelper.X1(t2)[1] - FunctionHelper.X1(t1)[1]) * FunctionHelper.VGamma1(t1)[1];
        double kernelDenominator = Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X1(t2)[1]), 2);
        
        return kernelNumerator * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t2)[0], 0, FunctionHelper.Der1X1(t2)[1], 0) /
               kernelDenominator;
    }

    double K12(double t1, double t2)
    {
        double kernelLeftValue = (FunctionHelper.VGamma1(t1)[0] * FunctionHelper.VGamma2(t2)[0] + FunctionHelper.VGamma1(t1)[1] * FunctionHelper.VGamma2(t2)[1]) /
                                 Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X2(t2)[1]), 2);
        double kernelRightValueNominator = ((FunctionHelper.X1(t1)[0] - FunctionHelper.X2(t2)[0]) * FunctionHelper.VGamma2(t2)[0] + (FunctionHelper.X1(t1)[1] - FunctionHelper.X2(t2)[1]) * FunctionHelper.VGamma2(t2)[1]) *
                                           ((FunctionHelper.X1(t1)[0] - FunctionHelper.X2(t2)[0]) * FunctionHelper.VGamma1(t1)[0] + (FunctionHelper.X1(t1)[1] - FunctionHelper.X2(t2)[1]) * FunctionHelper.VGamma1(t1)[1]);
        double kernelRightValue = kernelRightValueNominator / 
                                  Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X1(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X1(t1)[1], FunctionHelper.X2(t2)[1]), 4);
        
        return (kernelLeftValue - 2 * kernelRightValue) * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t2)[0], 0, FunctionHelper.Der1X2(t2)[1], 0);
    }

    double K21(double t1, double t2)
    {
        double kernelLeftValue = Math.Log(1 / FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X1(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X1(t2)[1]));
        double kernelRightValue = FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X1(t2)[0], 0, FunctionHelper.Der1X1(t2)[1], 0);
        return kernelLeftValue * kernelRightValue;
    }

    double K22(double t1, double t2)
    {
        if (Math.Abs(t1 - t2) < FunctionHelper.Eps)
        {
            return (-1 * FunctionHelper.Der1X2(t1)[0] * FunctionHelper.Der2X2(t1)[1] + FunctionHelper.Der1X2(t1)[1] * FunctionHelper.Der2X2(t1)[0]) /
                   (2 * Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t1)[0], 0, FunctionHelper.Der1X2(t1)[1], 0), 2));
        }

        double kernelNumerator = (FunctionHelper.X2(t1)[0] - FunctionHelper.X2(t2)[0]) * FunctionHelper.VGamma2(t2)[0] + (FunctionHelper.X2(t1)[1] - FunctionHelper.X2(t2)[1]) * FunctionHelper.VGamma2(t2)[1];
        double kernelDenominator = Math.Pow(FunctionHelper.GetEuclideanDistance(FunctionHelper.X2(t1)[0], FunctionHelper.X2(t2)[0], FunctionHelper.X2(t1)[1], FunctionHelper.X2(t2)[1]), 2);
        
        return kernelNumerator * FunctionHelper.GetEuclideanDistance(FunctionHelper.Der1X2(t2)[0], 0, FunctionHelper.Der1X2(t2)[1], 0) /
            kernelDenominator;
    }

    public double[] solutionValues = new double[2];
    
    public double[] Solve(int N1, double[] fValues, double[] gValues)
    {
        int N = N1;
        double[] H_F_Values = new double [4*N];
        double[,] kernelMatrix = new double [4*N, 4*N];
        for (var i = 0; i < N*4; i++) { H_F_Values[i] = i < 2*N ? fValues[i] : gValues[i - 2*N]; }
        // Console.WriteLine(H_F_Values[0] + " " +  H_F_Values[H_F_Values.Length - 1]);

        for (int i = 0; i < 2*N; i++)
        {
            double ti = i * Math.PI / N;
            for (int j = 0; j < 2*N; j++)
            {
                double tj = j * Math.PI / N;
                kernelMatrix[i, j] = K11(ti, tj)/(2*N);
                kernelMatrix[i, 2*N + j] = K12(ti, tj)/(2*N);
                kernelMatrix[2*N + i, j] = K21(ti, tj)/(2*N);
                kernelMatrix[2*N + i, 2*N + j] = K22(ti, tj)/(2*N);
            }

            kernelMatrix[i, i] -= 0.5;
            kernelMatrix[2*N + i, 2*N + i] -= 0.5;
        }

        double[,] kernelMatrixExtended = new double[4*N, 4*N + 1];
        for (int i = 0; i < 4*N; i++) { kernelMatrixExtended[i, 4*N] = H_F_Values[i]; }
        for (int i = 0; i < 4*N; i++) { for (int j = 0; j < 4*N; j++) { kernelMatrixExtended[i, j] = kernelMatrix[i, j]; } }
        
        // for (int i = 0; i < 4 * N; i++) { for (int j = 0; j < 4 * N + 1; j++) { Console.Write( kernelMatrixExtended[i,j]); } Console.WriteLine(); } Console.WriteLine();

        double[] ans = FunctionHelper.Gauss(kernelMatrixExtended, 4 * N);
        // Console.WriteLine("Gauss Values:");
        // for (int j = 0; j < 4 * N; j++) { Console.Write(ans[j] + " "); }
        // Console.WriteLine("\n");

        solutionValues = ans;
        return solutionValues;
    }
}