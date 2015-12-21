using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab5v2
{
    class Program
    {
        static void Main(string[] args)
        {
            int n = Convert.ToInt32(Console.ReadLine());
            Approximation App = new Approximation();
            Console.WriteLine("APolApp={0}", App.APolApprox(n));
            Console.WriteLine("ChPolApp={0}", App.ChPolApprox(n));
            Console.ReadKey();
        }
    }
}
