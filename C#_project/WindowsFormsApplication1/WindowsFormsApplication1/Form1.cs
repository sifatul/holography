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
using System.Drawing.Imaging;
using System.Diagnostics;

namespace WindowsFormsApplication1
{

    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private unsafe Bitmap ToBitmap(double[,] rawImage)
        {
            int width = rawImage.GetLength(1);
            int height = rawImage.GetLength(0);

            Bitmap Image = new Bitmap(width, height);
            BitmapData bitmapData = Image.LockBits(
                new Rectangle(0, 0, width, height),
                ImageLockMode.ReadWrite,
                PixelFormat.Format32bppArgb
            );
            ColorARGB* startingPosition = (ColorARGB*) bitmapData.Scan0;


            for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                {
                    double color = rawImage[i, j];
                    byte rgb = (byte)(color * 255);

                    ColorARGB* position = startingPosition + j + i * width;
                    position->A = 255;
                    position->R = rgb;
                    position->G = rgb;
                    position->B = rgb;
                }

            Image.UnlockBits(bitmapData);
            return Image;
        }

        public struct ColorARGB
        {
            public byte B;
            public byte G;
            public byte R;
            public byte A;

            public ColorARGB(Color color)
            {
                A = color.A;
                R = color.R;
                G = color.G;
                B = color.B;
            }

            public ColorARGB(byte a, byte r, byte g, byte b)
            {
                A = a;
                R = r;
                G = g;
                B = b;
            }

            public Color ToColor()
            {
                return Color.FromArgb(A, R, G, B);
            }
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
           

            int counter_y = 0;
            for (int i = 0; i < s; i++)
            {
                
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
            double[,] phase_h = new double[s, s];
            double max = -1.0;


            for (int i = 0; i < Cut.Count; i++)
            {
                film = new Complex[s, s];
                H = new Complex[s, s];
                O_image = new Complex[s, s];
                d1 = d - Cut[i] * Hologram_sampling_interval / 2;


                for (int j = 0; j < obj_z.Count; j++)
                {
                    if (Cut[i] == obj_z[j])
                    {
                        O_image[Convert.ToInt32(obj_x[j]), Convert.ToInt32(obj_y[j])] = (Complex)obj_r[j];
                    }
                }
                
                FourierTransform.FFT2(O_image, FourierTransform.Direction.Backward); // fft2    O_image = fft2(O_image); 
                
               // Stopwatch stopwatch = Stopwatch.StartNew(); //creates and start the instance of Stopwatch
               
                

                Complex com2 = Complex.Exp(new Complex(0, (1 * k * d1)));
                for (var p = 0; p < s; p++)
                {    //  H = exp(1i*k*d1).*exp(-1i*pi*lambda*d1*(x.^2+y.^2));

                    Complex com1;                    
                    for (var j = 0; j < s; j++)
                    {
                        com1 = Complex.Exp(new Complex(0, (-1 * Math.PI * lambda * d1 * (Math.Pow(x[p, j], 2) + Math.Pow(y[p, j], 2)))));
                       
                        H[p, j] = Complex.Multiply(com1, com2);  // H = exp(1i*k*d1).*exp(-1i*pi*lambda*d1*(x.^2+y.^2));
                        film[p, j] = Complex.Multiply(O_image[p, j], H[p, j]);  // lol =O_image.*H;
                    }

                }

                //////Stop();
                //Console.WriteLine(stopwatch.ElapsedMilliseconds);
                FourierTransform.FFT2(film, FourierTransform.Direction.Forward);    //ifft2    film3 =ifft2(lol); 

                
                for (int p = 0; p < s; p++)
                {

                    for (int j = 0; j < s; j++)
                    {
                        Hologram[p, j] = Complex.Add(Hologram[p, j], film[p, j]);  //  Hologram = Hologram+film3;

                        if (i == Cut.Count - 1) {

                            phase_h[p, j] = Hologram[p, j].Phase + Math.PI;    //phase_H = angle(Hologram) + pi;
                            if (phase_h[p, j] > max)
                            {
                                max = phase_h[p, j];
                            }
                        
                        
                        } // end of finding max and adding phase

                    }

                }



            } //end of Cut

         
               // Console.WriteLine(max);
            byte[,] phase_h_image = new byte[s,s];
            double d2 = d - o * 0.0001;
            
            Complex com3 = Complex.Exp(new Complex(0, (1 * k * -d2)));
            for (int p = 0; p <s; p++)
            {

                Complex com1;
               
                for (int j = 0; j < s; j++)
                {

                    double temp = 255 * phase_h[p, j] / max;   ///phase_H_image = uint8(255*phase_H/max(max(phase_H))); 
                    phase_h_image[p, j] = System.Convert.ToByte(temp); 
                    
                    com1 = Complex.Exp(new Complex(0, (-1 * Math.PI * lambda * -d2 * (Math.Pow(x[p, j], 2) + Math.Pow(y[p, j], 2)))));
                    H[p, j] = Complex.Multiply(com1, com3);         //H = exp(1i*k*z).*exp(-1i*pi*lambda*z*(x.^2+y.^2));

                }


            }
         
     
            FourierTransform.FFT2(Hologram, FourierTransform.Direction.Backward);  //  O = fft2(object);  

            Complex[,] originalR = new Complex[s, s];
            for (int p = 0; p < s; p++)
            {

                for (int j = 0; j < s; j++)
                {

                    originalR[p, j] = Complex.Multiply(Hologram[p, j], H[p, j]); //hologram =ifft2(O.*H); = >    O_F =  O.*H

                }

            }


            FourierTransform.FFT2(originalR, FourierTransform.Direction.Forward);  // hologram =ifft2(O_F)

            double[,] O_F = new double[s,s];
            for (int p = 0; p < s; p++)
             {

                 for (int j = 0; j <s; j++)
                 {
                     O_F[p, j] = originalR[p, j].Magnitude;  //abs(originalR)

                 }

             }


          Bitmap bmp = ToBitmap(O_F);
          PictureBox P = new PictureBox();


          P.Image = bmp;
          P.Dock = DockStyle.Fill;
          this.Controls.Add(P);
          this.Show();
  
       
            
            System.Console.WriteLine("END");
        }// end of form
        //System.Console.WriteLine("Finished");

    }
}// end of namespace

