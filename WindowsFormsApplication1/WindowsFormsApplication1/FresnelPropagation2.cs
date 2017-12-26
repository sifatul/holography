using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
namespace WindowsFormsApplication1
{
    class FresnelPropagation2
    {


        public FresnelPropagation2(Complex [][] hologram, double [,]x, double[,]y, double z, double lambda){




            Complex complex1 = new Complex(1,1);
            double k = 2 * Math.PI / lambda;
            Complex [][] H= new Complex[hologram.Length][];

            for (var p = 0; p < hologram.Length; p++)
                {    // fourier transform all values of x and y
                   H[0]= new Complex[hologram.Length];

                    Complex com1;
                    Complex com2;
                    for (var j = 0; j < hologram[p].Length; j++)
                    {
                        com1=  Complex.Exp( new Complex(0, (-complex1.Imaginary * Math.PI * lambda * z * (Math.Pow(x[p, j], 2) + Math.Pow(y[p, j], 2))) ) );
                        com2= Complex.Exp( new Complex(0, (complex1.Imaginary * k * z)));
                        H[p][j] = Complex.Multiply(com1, com2);


                       
                     }
                     
                }

          
        

        }


    }
}
