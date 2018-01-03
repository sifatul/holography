using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
namespace WindowsFormsApplication2
{
    class matrix_calculation_and_printing
    {



        public double[][] matrix_multiply(double[][] a, double [][] b) { 
        
            double [][] result= new double[a.Length][];

            for (int i = 0; i < b.Length; i++)
            {
                result[i] = new double[b.Length];

                for (int j = 0; j < b[0].Length; j++)
                {
                    double temp = 0;
                    for (int k = 0; k < a.Length; k++)
                    {

                        
                        temp += a[i][k] * b[k][j];

                    }

                    result[i][j]=temp;
                }
                //Console.WriteLine("\n");
            }

            return result;
            
        
        
        }



     
        public Complex[][] matrix_multiply(Complex[][] a, Complex[][] b)
        {

            Complex[][] result = new Complex[a.Length][];

            for (int i = 0; i < b.Length; i++)
            {
                result[i] = new Complex[b.Length];

                for (int j = 0; j < b[0].Length; j++)
                {
                    Complex temp = 0;
                    for (int k = 0; k < a.Length; k++)
                    {

                        
                        temp = Complex.Add(Complex.Multiply(a[i][k],b[k][j]),temp);

                    }

                    result[i][j] = temp;
                }
                //Console.WriteLine("\n");
            }

            return result;



        }


        public void print_matrix(double[][] b) {

            for (int i = 0; i < b.Length; i++)
            {
                

                for (int j = 0; j < b[0].Length; j++)
                {
                    Console.Write(b[i][j]+ " ");
                
                }
                Console.WriteLine();
            }
        
        
        }

        public Complex[][] add(Complex[][] a, Complex[][] b)
        {
            Complex[][] result = new Complex[a.Length][];

            for (int i = 0; i < b.Length; i++)
            {
                result[0] = new Complex[a.Length];

                for (int j = 0; j < b[0].Length; j++)
                {
                    result[i][j] = Complex.Add(a[i][j], b[i][j]);

                }
               
            }

            return result;

        }

        public void write_complex_in_file(String address, Complex[][] matrix) { 
        
            try
            {
                System.Console.WriteLine("Start writing");

                System.IO.StreamWriter file2 = new System.IO.StreamWriter(address, true);

                file2.WriteLine("              ################################                   ");
                for (var p = 0; p < matrix.Length; p++)
                {  
                    for (var j = 0; j < matrix[0].Length; j++)
                    {
                        file2.Write(matrix[p][j] + " ");
                            
                    }
                    file2.WriteLine();
                }
            }
            catch (Exception ex12)
            {
                System.Console.Write(ex12);
            }
        
        }



    }
}
