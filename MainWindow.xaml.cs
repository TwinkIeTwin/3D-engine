using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace _3D_engine
{
	public partial class MainWindow : Window
	{
		struct vect3D
		{
			public double x, y, z, w;
			public vect3D(double x, double y, double z)
			{
				this.x = x;
				this.y = y;
				this.z = z;
				this.w = z;
			}

			public vect3D(double x, double y, double z, double w)
			{
				this.x = x;
				this.y = y;
				this.z = z;
				this.w = w;
			}

			public static vect3D operator +(vect3D v1, vect3D v2)
			{
				return new vect3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
			}
			public static vect3D operator +(vect3D v1, double k)
			{
				return new vect3D(v1.x + k, v1.y + k, v1.z + k);
			}
			public static vect3D operator -(vect3D vect1, vect3D vect2)
			{
				return new vect3D(vect1.x - vect2.x, vect1.y - vect2.y, vect1.z - vect2.z);
			}
			public static vect3D operator -(vect3D vect1, double k)
			{
				return new vect3D(vect1.x - k, vect1.y - k, vect1.z - k);
			}
			public static vect3D operator *(vect3D vect1, vect3D vect2)
			{
				return new vect3D(vect1.x * vect2.x, vect1.y * vect2.y, vect1.z * vect2.z);
			}
			public static vect3D operator *(vect3D vect, double k)
			{
				return new vect3D(vect.x * k, vect.y * k, vect.z * k);
			}
			public static vect3D operator /(vect3D vect1, vect3D vect2)
			{
				return new vect3D(vect1.x / vect2.x, vect1.y / vect2.y, vect1.z / vect2.z);
			}
			public static vect3D operator /(vect3D vect, double k)	
			{
				return new vect3D(vect.x / k, vect.y / k, vect.z / k);
			}
			public double DotProduct(vect3D vect)
			{
				return this.x * vect.x + this.y * vect.y + this.z * vect.z;
			}
			public double GetLength()
			{
				return Math.Sqrt(this.DotProduct(this));
			}
			public void Normalize()
			{
				double l = this.GetLength();
				this.x /= l;
				this.y /= l;
				this.z /= l;
			}
			public vect3D CrossProduct(vect3D vect)
			{
				vect3D v = new vect3D();
				v.x = this.y * vect.z - this.z * vect.y;
				v.y = this.z * vect.x - this.x * vect.z;
				v.z = this.x * vect.y - this.y * vect.x;
				return v;
			}

		}
		SolidColorBrush getColor(double c)
		{
			byte a = Convert.ToByte(255 * Math.Abs(c));
			return new SolidColorBrush(System.Windows.Media.Color.FromRgb(a, a, a));
		}

		class triangle : IComparable<triangle>
		{
			public vect3D[] vect;

			public double luminosity;
			public triangle(vect3D v1, vect3D v2, vect3D v3)
			{
				vect = new vect3D[3];
				vect[0] = v1;
				vect[1] = v2;
				vect[2] = v3;
			}
			public triangle()
			{
				vect = new vect3D[3];
			}
			public int CompareTo(triangle other)
			{
				double thisAvvZ = (this.vect[0].z + this.vect[1].z + this.vect[2].z) / 3.0;
				double otherAvvZ = (other.vect[0].z + other.vect[1].z + other.vect[2].z) / 3.0;
				return otherAvvZ > thisAvvZ ? 1 : -1;
			}
		}

		matrix4x4 getMatrixRotationZ(double angle)
		{
			matrix4x4 matRotZ = new matrix4x4();
			matRotZ.m[0, 0] = Math.Cos(angle);
			matRotZ.m[0, 1] = Math.Sin(angle);
			matRotZ.m[1, 0] = -Math.Sin(angle);
			matRotZ.m[1, 1] = Math.Cos(angle);
			matRotZ.m[2, 2] = 1;
			matRotZ.m[3, 3] = 1;
			return matRotZ;
		}

		matrix4x4 getMatrixRotationX(double angle)
		{
			matrix4x4 matRotX = new matrix4x4();
			matRotX.m[0, 0] = 1;
			matRotX.m[1, 1] = Math.Cos(angle);
			matRotX.m[1, 2] = Math.Sin(angle);
			matRotX.m[2, 1] = -Math.Sin(angle);
			matRotX.m[2, 2] = Math.Cos(angle); ;
			matRotX.m[3, 3] = 1;
			return matRotX;
		}

		matrix4x4 getMatrixTranslation(double x, double y, double z)
		{
			matrix4x4 matrixTranslation = new matrix4x4();
			matrixTranslation.m[0, 0] = 1.0;
			matrixTranslation.m[1, 1] = 1.0;
			matrixTranslation.m[2, 2] = 1.0;
			matrixTranslation.m[3, 3] = 1.0;
			matrixTranslation.m[3, 0] = x;
			matrixTranslation.m[3, 1] = y;
			matrixTranslation.m[3, 2] = z;
			return matrixTranslation;
		}

		matrix4x4 getMatrixProjection(double fFovRad, double fAspectRation, double fNear, double fFar)
		{
			matrix4x4 matProj = new matrix4x4();
			matProj.m[0, 0] = fAspectRation * fFovRad;
			matProj.m[1, 1] = fFovRad;
			matProj.m[2, 2] = fFar / (fFar - fNear);
			matProj.m[3, 2] = (-fFar * fNear) / (fFar - fNear);
			matProj.m[2, 3] = 1.0;
			matProj.m[3, 3] = 0.0;
			return matProj;
		}

		class matrix4x4
		{
			public double[,] m;

			public static matrix4x4 operator *(matrix4x4 m1, matrix4x4 m2)
			{
				matrix4x4 res = new matrix4x4();
				for (int c = 0; c < 4; c++)
				{
					for (int r = 0; r < 4; r++)
					{
						res.m[r, c] = m1.m[r, 0] * m2.m[0, c] + m1.m[r, 1] * m2.m[1, c] + m1.m[r, 2] * m2.m[2, c] + m1.m[r, 3] * m2.m[3, c];
					}
				}
				return res;
			}
			public static vect3D operator *(matrix4x4 m, vect3D v)
			{
				vect3D res = new vect3D();
				res.x = v.x * m.m[0, 0] + v.y * m.m[1, 0] + v.z * m.m[2, 0] + v.w * m.m[3, 0];
				res.y = v.x * m.m[0, 1] + v.y * m.m[1, 1] + v.z * m.m[2, 1] + v.w * m.m[3, 1];
				res.z = v.x * m.m[0, 2] + v.y * m.m[1, 2] + v.z * m.m[2, 2] + v.w * m.m[3, 2];
				res.w = v.x * m.m[0, 3] + v.y * m.m[1, 3] + v.z * m.m[2, 3] + v.w * m.m[3, 3];
				return res;
			}
			public matrix4x4()
			{
				m = new double[4, 4];
			}
		}

		class mesh
		{
			public List<triangle> tris;

			public mesh()
			{
				tris = new List<triangle>();
			}

			// succesfull?
			public void loadFromObjFile(string fileName)
			{
				string[] file = System.IO.File.ReadAllLines(fileName);
				List<vect3D> verts = new List<vect3D>();

				// line = typeOfLine + v1 + v2 + v3
				foreach (string line in file)
				{
					string[] lineValues = line.Split(' ');
					char typeOfLine = line[0];
					//vertex
					if (typeOfLine == 'v')
					{
						vect3D vect = new vect3D(Convert.ToDouble(lineValues[1], System.Globalization.CultureInfo.InvariantCulture),
							Convert.ToDouble(lineValues[2], System.Globalization.CultureInfo.InvariantCulture),
							Convert.ToDouble(lineValues[3], System.Globalization.CultureInfo.InvariantCulture));
						verts.Add(vect);
					}
					//triangle
					else if (typeOfLine == 'f')
					{
						int v1, v2, v3;
						// - 1 because all of the information in the file start counting from 1
						v1 = Convert.ToInt32(lineValues[1]) - 1;
						v2 = Convert.ToInt32(lineValues[2]) - 1;
						v3 = Convert.ToInt32(lineValues[3]) - 1;
						tris.Add(new triangle(verts[v1], verts[v2], verts[v3]));
					}
				}
				
			}

		}
		
		matrix4x4 createUnitMatrix()
		{
			matrix4x4 m = new matrix4x4();
			m.m[0, 0] = 1.0;
			m.m[1, 1] = 1.0;
			m.m[2, 2] = 1.0;
			m.m[3, 3] = 1.0;
			return m;
		}
		void MultiplyMatrixVector(vect3D inputVect, out vect3D outputVect, matrix4x4 mat)
		{
			outputVect.x = inputVect.x * mat.m[0, 0] + inputVect.y * mat.m[1, 0] + inputVect.z * mat.m[2, 0] + mat.m[3, 0];
			outputVect.y = inputVect.x * mat.m[0, 1] + inputVect.y * mat.m[1, 1] + inputVect.z * mat.m[2, 1] + mat.m[3, 1];
			outputVect.z = inputVect.x * mat.m[0, 2] + inputVect.y * mat.m[1, 2] + inputVect.z * mat.m[2, 2] + mat.m[3, 2];
			outputVect.w = inputVect.x * mat.m[0, 3] + inputVect.y * mat.m[1, 3] + inputVect.z * mat.m[2, 3] + mat.m[3, 3];
			//deviding by zero
			if (outputVect.w != 0.0)
			{
				outputVect.x /= outputVect.w;
				outputVect.y /= outputVect.w;
				outputVect.z /= outputVect.w;
			}
		}

		void DrawTriangle(triangle t)
		{
			int x1 = Convert.ToInt32(t.vect[0].x);
			int y1 = Convert.ToInt32(t.vect[0].y);
			int x2 = Convert.ToInt32(t.vect[1].x);
			int y2 = Convert.ToInt32(t.vect[1].y);
			int x3 = Convert.ToInt32(t.vect[2].x);
			int y3 = Convert.ToInt32(t.vect[2].y);
			SolidColorBrush color = getColor(t.luminosity);

			System.Windows.Shapes.Polygon tri = new System.Windows.Shapes.Polygon();
			tri.Fill = color;
			tri.Points.Add(new System.Windows.Point(x1, y1));
			tri.Points.Add(new System.Windows.Point(x2, y2));
			tri.Points.Add(new System.Windows.Point(x3, y3));
			canvas.Children.Add(tri);
		}


		static double fTheta = 0;
		static double speed = 0.05;
		System.Windows.Threading.DispatcherTimer timer;

		private void initializeObjects()
		{
			// порядок указания вершин - по часовой стрелке
			meshCube = new mesh();
			meshCube.loadFromObjFile("spaceShip.obj");
		}

		private void initializeMatrix()
		{
			// min visability
			double fNear = 0.1;
			// max visible range
			double fFar = 1000;
			//field of view (in degrees)
			double fFov = 90;
			double fAspectRation = canvas.Height / canvas.Width; // Ww/Wh [readonly]
			// in rads
			double fFovRad = 1.0 / Math.Tan(fFov * Math.PI / 360);
			matProj = getMatrixProjection(fFovRad, fAspectRation, fNear, fFovRad);
		}

		public MainWindow()
		{
			vCamera = new vect3D(0, 0, 0);
			InitializeComponent();
			initializeMatrix();
			initializeObjects();
			System.Windows.Threading.DispatcherTimer timer = new System.Windows.Threading.DispatcherTimer();
			timer.Tick += new EventHandler(timerTick);
			timer.Interval = new TimeSpan(0, 0, 0, 0, 20);
			timer.Start();
		}
		mesh meshCube;
		matrix4x4 matProj;
		vect3D vCamera;
		private void timerTick(object sender, EventArgs e)
		{
			canvas.Children.Clear();

			fTheta += speed;
			matrix4x4 matRotZ = getMatrixRotationZ(fTheta);
			matrix4x4 matRotX = getMatrixRotationX(fTheta);

			matrix4x4 matTrans = getMatrixTranslation(0, 0, 16); ;

			// all the transformations that we need to do in this matrix
			matrix4x4 matWorld = createUnitMatrix();
			matWorld = matRotZ * matRotX;
			matWorld = matWorld * matTrans;

			List<triangle> trianglesToRaster = new List<triangle>();

			foreach (triangle tri in meshCube.tris)
			{
				triangle triProjected = new triangle();
				triangle triTransleted = new triangle();
				triangle triRotatedZ = new triangle();
				triangle triRotatedZX = new triangle();

				//test
				triangle triTransformed = new triangle();
				triTransformed.vect[0] = matWorld * tri.vect[0];
				triTransformed.vect[1] = matWorld * tri.vect[1];
				triTransformed.vect[2] = matWorld * tri.vect[2];


				//

				//rotate in Z
				MultiplyMatrixVector(tri.vect[0], out triRotatedZ.vect[0], matRotZ);
				MultiplyMatrixVector(tri.vect[1], out triRotatedZ.vect[1], matRotZ);
				MultiplyMatrixVector(tri.vect[2], out triRotatedZ.vect[2], matRotZ);

				// rotate in X
				MultiplyMatrixVector(triRotatedZ.vect[0], out triRotatedZX.vect[0], matRotX);
				MultiplyMatrixVector(triRotatedZ.vect[1], out triRotatedZX.vect[1], matRotX);
				MultiplyMatrixVector(triRotatedZ.vect[2], out triRotatedZX.vect[2], matRotX);

				// zoom out
				triTransleted = triRotatedZX;
				triTransleted.vect[0].z = triRotatedZX.vect[0].z + 8.0;
				triTransleted.vect[1].z = triRotatedZX.vect[1].z + 8.0;
				triTransleted.vect[2].z = triRotatedZX.vect[2].z + 8.0;

				//finding triangles normals
				vect3D normal = new vect3D();
				vect3D line1 = new vect3D();
				vect3D line2 = new vect3D();

				line1 = triTransleted.vect[1] - triTransleted.vect[0];
				line2 = triTransleted.vect[2] - triTransleted.vect[0];
				//line1.x = triTransleted.vect[1].x - triTransleted.vect[0].x;
				//line1.y = triTransleted.vect[1].y - triTransleted.vect[0].y;
				//line1.z = triTransleted.vect[1].z - triTransleted.vect[0].z;

				//line2.x = triTransleted.vect[2].x - triTransleted.vect[0].x;
				//line2.y = triTransleted.vect[2].y - triTransleted.vect[0].y;
				//line2.z = triTransleted.vect[2].z - triTransleted.vect[0].z;

				normal = line1.CrossProduct(line2);

				//normalizing normal
				normal.Normalize();


				//projecting

				// объект виден если угол между нормалью к нему и нашей камерой меньше 90 градусов
				// будем проверять это с помощью скалярного произведения векторов
				// скалярное произведение векторов угол между которыми > 90 градусов будет < 0
				if (normal.DotProduct(triTransleted.vect[0] - vCamera) < 0.0)
				{
					// illumination (all illumination cooming from direction, not from the point)
					// свет тем ярче чем меньше угол между нормалью поверхности с направлением света
					vect3D lightDirection = new vect3D(0.0, 0.0, -1.0);
					lightDirection.Normalize();

					// dot product between normal and lightDirection
					double dp = normal.x * lightDirection.x + normal.y * lightDirection.y + normal.z * lightDirection.z;

					triProjected.luminosity = dp;


					MultiplyMatrixVector(triTransleted.vect[0], out triProjected.vect[0], matProj);
					MultiplyMatrixVector(triTransleted.vect[1], out triProjected.vect[1], matProj);
					MultiplyMatrixVector(triTransleted.vect[2], out triProjected.vect[2], matProj);



					//Scale into view
					triProjected.vect[0].x += 1.0;
					triProjected.vect[0].x *= SystemParameters.PrimaryScreenWidth / 2.0;
					triProjected.vect[0].y += 1.0;
					triProjected.vect[0].y *= SystemParameters.PrimaryScreenHeight / 2.0;

					triProjected.vect[1].x += 1.0;
					triProjected.vect[1].x *= SystemParameters.PrimaryScreenWidth / 2.0;
					triProjected.vect[1].y += 1.0;
					triProjected.vect[1].y *= SystemParameters.PrimaryScreenHeight / 2.0;

					triProjected.vect[2].x += 1.0;
					triProjected.vect[2].x *= SystemParameters.PrimaryScreenWidth / 2.0;
					triProjected.vect[2].y += 1.0;
					triProjected.vect[2].y *= SystemParameters.PrimaryScreenHeight / 2.0;

					trianglesToRaster.Add(triProjected);
				}
			}
			// Sort triangles from back to front 
			//MessageBox.Show(trianglesToRaster[0].vect[0].x + " " + trianglesToRaster[0].vect[0].y + " " + trianglesToRaster[0].vect[0].z);
			trianglesToRaster.Sort();
			//MessageBox.Show(trianglesToRaster[0].vect[0].x + " " + trianglesToRaster[0].vect[0].y + " " + trianglesToRaster[0].vect[0].z);
			foreach (triangle t in trianglesToRaster)
			{
				DrawTriangle(t);
			}
		}

	}
}
