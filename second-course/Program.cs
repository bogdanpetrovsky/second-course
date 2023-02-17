﻿// See https://aka.ms/new-console-template for more information
double[] X1(double t) { return new [] { Math.Cos(t), Math.Sin(t) }; }
double[] X2(double t) { return new [] { 2 * Math.Cos(t), 2 * Math.Sin(t) }; }
double[] Der1X1(double t) { return new [] { -1 * Math.Sin(t), Math.Cos(t) }; }
double[] Der1X2(double t) { return new [] { -2 * Math.Sin(t), 2 * Math.Cos(t) }; }
double[] Der2X1(double t) { return new [] { -1 * Math.Cos(t), -1 * Math.Sin(t) }; }
double[] Der2X2(double t) { return new [] { -2 * Math.Cos(t), -2 * Math.Sin(t) }; }

double[] Vy(double t)
{
    return new []
    {
        Der1X2(t)[1] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0),
        -1 * Der1X2(t)[0] / GetEuclideanDistance(Der1X2(t)[0], 0, Der1X2(t)[1], 0)
    };
}
double[] Vx(double t)
{
    return new []
    {
        Der1X1(t)[1] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0),
        -1 * Der1X1(t)[0] / GetEuclideanDistance(Der1X1(t)[0], 0, Der1X1(t)[1], 0)
    };
}

double Kernel11(double t1, double t2)
{
    if (t1 == t2)
    {
        return (-1 * Der1X1(t1)[0] * Der2X1(t1)[1] + Der1X1(t1)[1] * Der2X1(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0), 2));
    }

    double kernelNumerator = (X1(t1)[0] - X1(t2)[0]) * Vx(t1)[0] + (X1(t1)[1] - X1(t2)[1]) * Vx(t1)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X1(t1)[0], X1(t2)[0], X1(t1)[1], X1(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X1(t1)[0], 0, Der1X1(t1)[1], 0) /
           kernelDenominator;
}

double Kernel12(double t1, double t2)
{
    double kernelLeftValue = (Vy(t2)[0] * Vx(t1)[0] + Vy(t2)[1] * Vx(t1)[1]) /
                             GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]);
    double kernelRightValueNominator = (X1(t1)[0] - X2(t2)[0] * Vy(t2)[0] + X1(t1)[1] - X2(t2)[1] * Vy(t2)[1]) *
                                       (X1(t1)[0] - X2(t2)[0] * Vx(t1)[0] + X1(t1)[1] - X2(t2)[1] * Vx(t1)[1]);
    double kernelRightValue = kernelRightValueNominator / 
                              Math.Pow(GetEuclideanDistance(X1(t1)[0], X2(t2)[0], X1(t1)[1], X2(t2)[1]), 4);
    
    return 0.5 * Math.PI * (kernelLeftValue - 2 * kernelRightValue) * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0);
}

double Kernel21(double t1, double t2)
{
    double kernelLeftValue = Math.Log(1 / GetEuclideanDistance(X2(t1)[0], X1(t2)[0], X2(t1)[1], X1(t2)[1]));
    double kernelRightValue = GetEuclideanDistance(Der1X1(t2)[0], 0, Der1X1(t2)[1], 0);
    return kernelLeftValue * kernelRightValue;
}

double Kernel22(double t1, double t2)
{
    if (t1 == t2)
    {
        return (-1 * Der1X2(t1)[0] * Der2X2(t1)[1] + Der1X2(t1)[1] * Der2X2(t1)[0]) /
               (2 * Math.Pow(GetEuclideanDistance(Der1X2(t1)[0], 0, Der1X2(t1)[1], 0), 2));
    }

    double kernelNumerator = (X2(t1)[0] - X2(t2)[0]) * Vy(t2)[0] + (X2(t1)[1] - X2(t2)[1]) * Vy(t2)[1];
    double kernelDenominator = Math.Pow(GetEuclideanDistance(X2(t1)[0], X2(t2)[0], X2(t1)[1], X2(t2)[1]), 2);
    
    return kernelNumerator * GetEuclideanDistance(Der1X2(t2)[0], 0, Der1X2(t2)[1], 0) /
        kernelDenominator;
}

double GetEuclideanDistance(double x21, double x11, double x22, double x12)
{
    return Math.Sqrt(Math.Pow(x21 - x11, 2) + Math.Pow(x22 - x12, 2));
}

void Solve()
{
    int N = 4;
    double[] H_F_Values = new double [4*N];
    double[,] kernelMatrix = new double [4*N, 4*N];
    for (var i = 0; i < N * 4; i++) { H_F_Values[i] = i < 2 * N ? 0 : 1; }

    for (int i = 0; i < 4*N; i++)
        for (int j = 0; j < 4*N; j++)
        {
            if (i == j && i < 2*N) { kernelMatrix[i, j] = 0.5 * N * Kernel11(i * Math.PI/N, j * Math.PI/N) - 0.5 * i; }
            if (i == j && i >= 2*N) { kernelMatrix[i, j] = 0.5 * N * Kernel11(i * Math.PI/N, j * Math.PI/N) - 0.5 * i; }
            if (i < 2*N && j < 2*N) { kernelMatrix[i, j] = Kernel11(i * Math.PI/N, j * Math.PI/N); }
            if (i < 2*N && j >= 2*N) { kernelMatrix[i, j] = Kernel12(i * Math.PI/N, j * Math.PI/N); }
            if (i >= 2*N && j < 2*N) { kernelMatrix[i, j] = Kernel21(i * Math.PI/N, j * Math.PI/N); }
            if (i >= 2*N && j >= 2*N) { kernelMatrix[i, j] = Kernel22(i * Math.PI/N, j * Math.PI/N); }
        }
    
    for (int i = 0; i < 4 * N; i++)
    {
        for (int j = 0; j < 4 * N; j++)
        {
            Console.Write(kernelMatrix[i,j] + " ");
        }
        Console.WriteLine();
    }
    Console.WriteLine();
    for (int j = 0; j < 4 * N; j++)
    {
        Console.Write(H_F_Values[j] + " ");
    }
    
    // call Gauss
}

Solve();
Console.WriteLine("\nEND");