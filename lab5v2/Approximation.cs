using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace lab5v2
{
    class Approximation
    {
        interface IIntegrate
        {
            double Func(double x);
        }
        double[] MainElement(int N, double[,] M)
        {
            int[] position = new int[N];
            for (int i = 0; i < N; i++)
                position[i] = i;
            for (int k =N;k>1;k--)
            {
                double max = M[0, 0];
                double[] m = new double[k];
                
                int r = 0, c = 0;
                for (int i=0;i< k;i++)
                {
                    for(int j=0;j< k;j++)
                        if(M[i,j]>max)
                        {
                            max = M[i, j];
                            r = i;
                            c = j;
                        }
                }
                for(int i = 0; i< k;i++)
                {
                    m[i] = -M[i, c] / max;
                    if (i == r)
                        m[i] = 0;
                }
                for(int i=0;i< k;i++)
                {
                    for(int j=0;j< k;j++)
                    {
                        M[i, j] += M[r, j] * m[i];
                    }
                    M[i, N] += M[r, N] * m[i];
                }
                for(int i=0;i<N+1;i++)
                {
                    double d1 = M[r, i];
                    M[r, i] = M[k - 1, i];
                    M[k - 1, i] = d1;
                }
                for(int i=0;i< N;i++)
                {
                    double d1 = M[i, c];
                    M[i, c] = M[i, k - 1];
                    M[i, k - 1] = d1;
                }
                int d = position[c];
                position[c] = position[k - 1];
                position[k - 1] = d;
            }
            double[] x1 = new double[N];
            for(int i=0;i< N;i++)
            {
                x1[i] = 0;
            }
            double[] x2 = new double[N];
            for(int i=0;i<N;i++)
            {
                double sum = M[i, N];
                for(int j=0;j< i;j++)
                {
                    sum -= M[i, j] * x1[j];
                }
                x1[i] = sum / M[i, i];
            }
            for (int i = 0; i < N; i++)
                x2[position[i]] = x1[i];

            return x2;
        }
        double Integr(IIntegrate F,int n,double a, double b)
        {
            
            double h = (b - a) / n;
            double y0 = F.Func(a);
            double y1 = F.Func(b);
            double s1 = 0, s2 = 0;
            double x = a + h;
            for(int i=1;i< n;i++)
            {
                double y = F.Func(x);
                if (i % 2 == 1)
                    s1 += y;
                else
                    s2 += y;
                x += h;
                
            }
            return h * (y0 + y1 + 4 * s1 + 2 * s2) / 3.0;
        }
        double Integr_Rynge(IIntegrate F, double a, double b)
        {
            double Eps = 10e-7;
            int n = (int)Math.Round(1 / Math.Sqrt(Math.Sqrt(Eps)));
            double i1 = Integr(F, n, a, b);
            n *= 2;
            double i2 = Integr(F, n, a, b);
            while((Math.Abs((i1- i2)/(15))>Eps*i2)&&(i2>Eps))
            {
                i1 = i2;
                n *= 2;
                i2 = Integr(F, n, a, b);
            }
            return i2;
        }
        class APol:IIntegrate
        {
            int N;
            public APol(int n)
            {
                N = n;
            }
            public double Func(double x)
            {
                return Math.Pow(x, N);
            }
        }
        class APolFunc:IIntegrate
        {
            int N; 
            public APolFunc(int n)
            {
                N = n;
            }
            public double F(double x)
            {
                return Math.Log10(x * x) * Math.Sin(x / 2) * Math.Exp(Math.Pow(x, 1.0 / 7.0));
                //return Math.Abs(x);
            }
            public double Func(double x)
            {
                return F(x) * Math.Pow(x, N);
            }
        }

        class APolAv : IIntegrate
        {
            int N;
            double[] A;
            public APolAv(int n,  double[] a )
            {
                N = n;
                A = a;
            }
            public double F(double x)
            {
                return Math.Log10(x * x) * Math.Sin(x / 2) * Math.Exp(Math.Pow(x, 1.0 / 7.0));
                //return Math.Abs(x);
            }
            public double APol(double x)
            {
                double sum = 0;
                for (int i = 0; i < N; i++)
                    sum += A[i] * Math.Pow(x, i);
                return sum;
            }
            public double Func(double x)
            {
                return Math.Pow(F(x) - APol(x), 2);
            }
        }
        class ChPol2 : IIntegrate
        {
            int N;
            public ChPol2(int n)
            {
                N = n;
            }
            double Pol(int n, double x)
            {
                if (n == 0)
                    return 1;
                if (n == 1)
                    return x;
                return 2 * x * Pol(n - 1, x) - Pol(n - 2, x);
            }
            public double Func(double x)
            {
                return Math.Pow(Pol(N, x),2);
            }
        }
        class ChPolF : IIntegrate
        {
            int N;
            public ChPolF(int n)
            {
                N = n;
            }
            double Pol(int n, double x)
            {
                if (n == 0)
                    return 1;
                if (n == 1)
                    return x;
                return 2 * x * Pol(n - 1, x) - Pol(n - 2, x);
            }
            double F(double t)
            {
                double x = 5.5 + 4.5 * t;
                return Math.Log10(x * x) * Math.Sin(x / 2) * Math.Exp(Math.Pow(x, 1.0 / 7.0));
                //return Math.Abs(x);
            }
            public double Func(double x)
            {
                return F(x) * Pol(N, x);
            }

        }
        class ChPolFAv:IIntegrate
        {
            int N;
            double[] A;
            public ChPolFAv(int n, double[] a)
            {
                A = a;
                N = n;
            }
            double Pol(int n, double x)
            {
                if (n == 0)
                    return 1;
                if (n == 1)
                    return x;
                return 2 * x * Pol(n - 1, x) - Pol(n - 2, x);
            }
            double F(double t)
            {
                double x = 5.5 + 4.5 * t;
                return Math.Log10(x * x) * Math.Sin(x / 2) * Math.Exp(Math.Pow(x, 1.0 / 7.0));
                //return Math.Abs(x);
            }
            public double FPol(double x)
            {
                double sum = 0;
                for(int i=0;i< N;i++)
                {
                    sum += A[i] * Pol(i, x);
                }
                return sum;
            }
            public double Func(double x)
            {
                return Math.Pow(F(x) - FPol(x), 2);
            }
        }

        public double APolApprox(int N)
        {
            double a = 1, b = 10;
            double[,] M= new double[N, N + 1];
            for(int i=0;i< N;i++)
            {
                for(int j=0;j< N;j++)
                {
                    APol F = new APol(i + j);
                    M[i, j] = Integr_Rynge(F, a, b);
                }
                APolFunc AF = new APolFunc(i);
                M[i, N] = Integr_Rynge(AF, a, b);
                //M[i, N] = AF.Func(35);
            }
            double[] A = MainElement(N, M);
            APolAv PolAv = new APolAv(N, A);
            double h = (b - a) / 200;
            double x = a;
            using (StreamWriter SW = new StreamWriter("out.txt"))
            {
                SW.WriteLine("{0};{1}", x, PolAv.APol(x));
                while (x <= b)
                {
                    x += h;
                    SW.WriteLine("{0};{1}", x, PolAv.APol(x));

                }
            }
            return Math.Sqrt(Integr_Rynge(PolAv, a, b) / (b - a));
        }
        public double ChPolApprox(int N)
        {
            double[] A = new double[N];
            double a = -1, b = 1;
            for (int i = 0; i < N; i++)
            {
                ChPol2 CP = new ChPol2(i);
                ChPolF CPF = new ChPolF(i);
                A[i] = Integr_Rynge(CPF, a, b) / Integr_Rynge(CP, a, b);

            }

            ChPolFAv CPA = new ChPolFAv(N, A);
            double h = (b - a) / 200;
            double x = a;
            using (StreamWriter SW = new StreamWriter("out1.txt"))
            {
                SW.WriteLine("{0};{1}", 5.5+4.5*x, CPA.FPol(x));
                while (x <= b)
                {
                    x += h;
                    SW.WriteLine("{0};{1}", 5.5 + 4.5 * x, CPA.FPol(x));

                }
            }
            return Math.Sqrt(Integr_Rynge(CPA, a, b) / (b - a));
        }
    }
}
