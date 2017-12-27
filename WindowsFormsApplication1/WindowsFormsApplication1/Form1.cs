using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
//using System.Numerics;
using AForge.Math;
using System.IO;


namespace WindowsFormsApplication1
{






    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }



        private void Form1_Load(object sender, EventArgs e)
        {


            String address = "D:/workshop/holography/WindowsFormsApplication1/WindowsFormsApplication1";
            
            String file_name = "/Tang25.txt";

            file_manager file = new file_manager(file_name, address);
            file.read_coordinate();
            // file.write_coordinate_file(file.x, file.y, file.z, file.r, file.g, file.b);


            List<double> obj_x = new List<double>();
            List<double> obj_y = new List<double>();
            List<double> obj_z = new List<double>();
            List<double> obj_r = new List<double>();
            List<double> obj_g = new List<double>();
            List<double> obj_b = new List<double>();

            List<double> O = new List<double>();
            for (int i = 0; i < file.x.Count; i++)
            {

                //System.Console.WriteLine(i + " " + file.x[i] + " " + file.y[i] + " " + file.z[i]);

                obj_x.Add(Math.Round((Convert.ToDouble(file.x[i]) + 100) * 5));
                obj_y.Add(Math.Round((Convert.ToDouble(file.y[i]) + 50) * 6));
                obj_z.Add(-Math.Round(Convert.ToDouble(file.z[i]) - 230));

                obj_r.Add(Convert.ToDouble(file.r[i]));
                obj_g.Add(Convert.ToDouble(file.g[i]));
                obj_b.Add(Convert.ToDouble(file.b[i]));


            }

            int o = 8;
            int t = 1;
            int s = 1024;
            double d = 0.25;
            double lambda = 532e-9;
            double k = 2 * Math.PI / lambda;

            double Hologram_sampling_interval = 7.4e-6;
            double dx = Hologram_sampling_interval;
            double dy = Hologram_sampling_interval;

            double[,] image = new double[s, s];



            for (int i = 0; i < obj_x.Count; i++)
            {
                obj_x[i] = obj_x[i] * t;
                obj_y[i] = obj_y[i] * t;
                obj_z[i] = obj_z[i] * 1;
            }


            List<double> Cut = obj_z.Distinct().ToList();
            Cut.Sort();
            double Ny = s;
            double Nx = s;
            double fx = 1 / (dx * Nx);
            double fy = 1 / (dy * Ny);

            double[,] x = new double[s, s];
            double[,] y = new double[s, s];





            Complex complex1 = new Complex(1, 1);

           
            
            Complex[,] Hologram = new Complex[s,s];
            Complex[][] O_F;

            int counter_y = 0;
            for (int i = 0; i < s; i++)
            {
                // O_image[i]=new Complex[s];
                // film[i] = new Complex[s];
                //H[i] = new Complex[s];
              //  Hologram[i] = new Complex[s];



                int counter_x = 0;
                if (counter_y == s / 2)
                {
                    counter_y *= -1;
                }

                for (int j = 0; j < s; j++)
                {

                    if (counter_x == s / 2)
                    {
                        counter_x *= -1;
                    }

                    x[i, j] = counter_x * fx;
                    y[i, j] = counter_y * fy;

                    counter_x++;

                }

                counter_y++;

            }


            double d1;

            Complex[,] O_image;
            Complex[,] H = new Complex[s, s]; 
            Complex[,] film;

            for (int i = 0; i < Cut.Count; i++)
            {
                film = new Complex[s, s];
                H = new Complex[s, s];
                O_image = new Complex[s, s];


                for (int j = 0; j < obj_z.Count; j++)
                {
                    if (Cut[i] == obj_z[j])
                    {

                        /* double[] temp = new double[6];
                         temp[0] = obj_x[j];
                         temp[1] = obj_y[j];
                         temp[2] = obj_z[j];

                         temp[3] = obj_r[j];
                         temp[4] = obj_g[j];
                         temp[5] = obj_b[j];
                         O.AddRange(temp);
                         */
                        O_image[Convert.ToInt32(obj_x[j]), Convert.ToInt32(obj_y[j])] = (Complex)obj_r[j];
                    }
                }

                d1 = d - Cut[i] * Hologram_sampling_interval / 2;


                FourierTransform.FFT2(O_image, FourierTransform.Direction.Backward); // fft2


                for (var p = 0; p < s; p++)
                {    // fourier transform all values of x and y

                    Complex com1;
                    Complex com2;
                    for (var j = 0; j < s; j++)
                    {
                        com1 = Complex.Exp(new Complex(0, (-complex1.Im * Math.PI * lambda * d1 * (Math.Pow(x[p, j], 2) + Math.Pow(y[p, j], 2)))));
                        com2 = Complex.Exp(new Complex(0, (complex1.Im * k * d1)));
                        H[i, j] = Complex.Multiply(com1, com2);

                    }

                }

                for (int p = 0; p < s; p++)
                {

                    for (int j = 0; j < s; j++)
                    {

                        film[i, j] = Complex.Multiply(O_image[i, j], H[i, j]);

                    }

                }

                FourierTransform.FFT2(film, FourierTransform.Direction.Forward);    //ifft2


                //until now the result is matched

                for (int p = 0; p < s; p++)
                {


                    for (int j = 0; j < s; j++)
                    {
                        Hologram[p, j] = Complex.Add(Hologram[p, j], film[p, j]);

                    }

                }



            } //end of Cut

            //hologram is not accurate match
            double[,] phase_h = new double[s, s];
            double max = -1.0;
            for (int p = 0; p < s; p++)
            {
                for (int j = 0; j <s; j++)
                {

                    phase_h[p, j] = Hologram[p,j].Phase + Math.PI;
                    if (phase_h[p, j] > max)
                    {
                        max = phase_h[p, j];
                    }
                }
            }

           // Console.WriteLine(max);
            double[,] phase_h_image = new double[s,s];
            for (int p = 0; p <s; p++)
            {
               

                for (int j = 0; j < s; j++)
                {

                    double temp = 255 * phase_h[p, j] / max;
                    phase_h_image[p,j] = System.Convert.ToByte(temp);

                }


            }

            double d2 = d - o * 0.0001;


            /*FresnelPropagation2

             * 
             * 
             * 
             * */



            for (var p = 0; p < s; p++)
            {    // fourier transform all values of x and y

                Complex com1;
                Complex com2;
                for (var j = 0; j < s; j++)
                {
                    com1 = Complex.Exp(new Complex(0, (-complex1.Im * Math.PI * lambda * -d2 * (Math.Pow(x[p, j], 2) + Math.Pow(y[p, j], 2)))));
                    com2 = Complex.Exp(new Complex(0, (complex1.Im * k * -d2)));
                    H[p, j] = Complex.Multiply(com1, com2);

                }

            }


            // O_F = FourierTransform2.FFT2(Hologram, "Forward");

            /* for (int p = 0; p < O_F.Length; p++)
             {

                 for (int j = 0; j < O_F[0].Length; j++)
                 {
                     Hologram[p][j] = Complex.Multiply(O_F[p][j], H[p][j]);

                 }

             }
             * /
             //Hologram = FourierTransform2.FFT2(film, "Backward");




            
             this.BackColor = Color.FromArgb(255, 102, 178);
         */

            System.Console.WriteLine("END");
        }// end of form
        //System.Console.WriteLine("Finished");

    }
}// end of namespace


/*

                file_name = "fft.txt";
                try
                {
                    System.Console.WriteLine("Start writing");


                    System.IO.StreamWriter file2 = new System.IO.StreamWriter(address + file_name, true);
                    for (var p = 0; p < s; p++)
                    {    // fourier transform all values of x and y

                        for (var j = 0; j < s; j++)
                        {


                            file2.Write(O_F[p][j] + " ");
                           

                            // file2.WriteLine(" ########################################################################### ");
                        }
                    }
                }
                catch (Exception ex12)
                {
                    System.Console.Write(ex12);
  */
