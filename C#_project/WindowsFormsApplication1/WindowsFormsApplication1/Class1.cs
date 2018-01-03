using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WindowsFormsApplication1
{
    public class file_manager
    {
        String file_name;
        String file_address;


        public List<string> x ;
        public List<string> y;
        public  List<string> z ;
        public List<string> r ;
        public List<string> g;
        public List<string> b ;

        public file_manager(String name, String address)
        {
            file_name = name;
            file_address = address;
           // System.Console.WriteLine(name);




        }

        public void convert_and_write_coordinate_file(List< Double> x, List<Double> y, List<Double> z)
        {


         
      
        }

        public void write_coordinate_file(List<string> x, List<string> y, List<string> z, List<string> r, List<string> g, List<string> b)
        {     

            

             
            System.IO.File.WriteAllLines(file_address+  "x.txt", x);
            System.IO.File.WriteAllLines(file_address + "y.txt", y);
            System.IO.File.WriteAllLines(file_address + "z.txt", z);
            System.IO.File.WriteAllLines(file_address + "r.txt", r);
            System.IO.File.WriteAllLines(file_address + "g.txt", g);
            System.IO.File.WriteAllLines(file_address + "b.txt", b);
                

        }

        public void write_coordinate_file(List<string> x, List<string> y, List<string> z)
        {

            System.IO.File.WriteAllLines(file_address + "x.txt", x);
            System.IO.File.WriteAllLines(file_address + "y.txt", y);
            System.IO.File.WriteAllLines(file_address + "z.txt", z);
        }

        public void read_coordinate()
        {
           
           

           x = new List<string>();
           y = new List<string>();
           z = new List<string>();
           r = new List<string>();
           g = new List<string>();
           b = new List<string>();

            try
            {
                System.IO.StreamReader file= new System.IO.StreamReader(file_address + file_name);
                

                String line;
                while ((line = file.ReadLine()) != null)
                {


                    //System.Console.WriteLine("read frin fuke");

                    var temp = line.Split(' ');

                    x.Add(temp[1]);
                    y.Add(temp[2]);
                    z.Add(temp[3]);

                    r.Add(temp[4]);
                    g.Add(temp[5]);
                    b.Add(temp[6]);


                }

                file.Close();

            }
            catch (Exception ex)
            {
                System.Console.WriteLine(ex);
            }







        }



        
            
            



    }
}
