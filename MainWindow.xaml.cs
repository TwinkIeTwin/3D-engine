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
				this.w = 1;
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
		//FUNC
		SolidColorBrush getColor(double c)
		{
			byte a = Convert.ToByte(255 * Math.Abs(c));
			a = a < 25 ? (byte)25 : a;
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
		//FUNC
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
		matrix4x4 getMatrixRotationY(double angle)
		{
			matrix4x4 matrix = new matrix4x4();
			matrix.m[0,0] = Math.Cos(angle);
			matrix.m[0,2] = Math.Sin(angle);
			matrix.m[2,0] = -Math.Sin(angle);
			matrix.m[1,1] = 1.0f;
			matrix.m[2,2] = Math.Cos(angle);
			matrix.m[3,3] = 1.0f;
			return matrix;
		}
		//FUNC
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
		//FUNC
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
		//FUNC
		matrix4x4 matrixPointAt(vect3D pos, vect3D target, vect3D up)
		{
			vect3D newForward = target - pos;
			newForward.Normalize();
			// Calculate new Up direction
			vect3D a = newForward * up.DotProduct(newForward);
			vect3D newUp = up - a;
			newUp.Normalize();
			// New Right direction is easy, its just cross product
			vect3D newRight = newUp.CrossProduct(newForward);
			// Construct Dimensioning and Translation Matrix	
			matrix4x4 matrix = new matrix4x4();
			matrix.m[0, 0] = newRight.x; matrix.m[0,1] = newRight.y; matrix.m[0, 2] = newRight.z; matrix.m[0,3] = 0.0f;
			matrix.m[1, 0] = newUp.x; matrix.m[1,1] = newUp.y; matrix.m[1,2] = newUp.z; matrix.m[1,3] = 0.0f;
			matrix.m[2, 0] = newForward.x; matrix.m[2,1] = newForward.y; matrix.m[2,2] = newForward.z; matrix.m[2,3] = 0.0f;
			matrix.m[3, 0] = pos.x; matrix.m[3,1] = pos.y; matrix.m[3,2] = pos.z; matrix.m[3,3] = 1.0f;
			return matrix;
		}

		matrix4x4 matrixQuickInverse(matrix4x4 m)
		{
			matrix4x4 matrix = new matrix4x4();
			matrix.m[0, 0] = m.m[0, 0];
			matrix.m[0, 1] = m.m[1, 0];
			matrix.m[0, 2] = m.m[2, 0];
			matrix.m[0, 3] = 0.0f;

			matrix.m[1, 0] = m.m[0, 1];
			matrix.m[1, 1] = m.m[1, 1];
			matrix.m[1, 2] = m.m[2, 1];
			matrix.m[1, 3] = 0.0f;

			matrix.m[2, 0] = m.m[0, 2];
			matrix.m[2, 1] = m.m[1, 2];
			matrix.m[2, 2] = m.m[2, 2];
			matrix.m[2, 3] = 0.0f;

			matrix.m[3, 0] = -(m.m[3, 0] * matrix.m[0, 0] + m.m[3, 1] * matrix.m[1, 0] + m.m[3, 2] * matrix.m[2, 0]);
			matrix.m[3, 1] = -(m.m[3, 0] * matrix.m[0, 1] + m.m[3, 1] * matrix.m[1, 1] + m.m[3, 2] * matrix.m[2, 1]);
			matrix.m[3, 2] = -(m.m[3, 0] * matrix.m[0, 2] + m.m[3, 1] * matrix.m[1, 2] + m.m[3, 2] * matrix.m[2, 2]);
			matrix.m[3, 3] = 1.0f;
			return matrix;
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
			//FUNC
			public void loadFromObjFileAt(string fileName, double x, double y, double z)
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
						vect3D vect = new vect3D(Convert.ToDouble(lineValues[1], System.Globalization.CultureInfo.InvariantCulture) + x,
							Convert.ToDouble(lineValues[2], System.Globalization.CultureInfo.InvariantCulture) + y,
							Convert.ToDouble(lineValues[3], System.Globalization.CultureInfo.InvariantCulture) + z);
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
		//FUNC
		void DrawTriangle(triangle t, bool points, bool lines)
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
			if (lines)
			{
				tri.Stroke = System.Windows.Media.Brushes.Red;
			}
			else
			{
				tri.Stroke = color;
			}
			tri.Points.Add(new System.Windows.Point(x1, y1));
			tri.Points.Add(new System.Windows.Point(x2, y2));
			tri.Points.Add(new System.Windows.Point(x3, y3));
			canvas.Children.Add(tri);
			if (points)
			{
				foreach (var i in t.vect)
				{
					Ellipse circle = new Ellipse();
					circle.Margin = new Thickness(i.x, i.y, 0, 0);
					circle.Width = 5;
					circle.Height = 5;
					circle.Fill = System.Windows.Media.Brushes.Blue;
					if (!canvas.Children.Contains(circle))
					{
						canvas.Children.Add(circle);
					}
				}
			}

		}

		static double fTheta = Math.PI / 2.0;
		static double speed = 0.05;
		System.Windows.Threading.DispatcherTimer timer;

		private void initializeObjects()
		{
			meshCube = new mesh();
			meshCube.loadFromObjFileAt("spaceShip.obj", 0, 0, 0);
			meshCube.loadFromObjFileAt("spaceShip.obj", 0, 2, 15);
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
			timer.Interval = new TimeSpan(0, 0, 0, 0, 1);
			timer.Start();
			canvas.Focus();
		}
		mesh meshCube;
		matrix4x4 matProj;
		vect3D vCamera;
		vect3D vLookDir;
		// camera rotation
		// Y axis (horizontal)
		double fYaw;
		// X axis (vertical)
		double fXaw;

		double x = 0;
		double r = 9;

		bool isHalfRotated = false;
		private void timerTick(object sender, EventArgs e)
		{
			canvas.Children.Clear();
	
			double deltaX = 0.2 ;
			double fDelta = 2 * Math.PI / (4 * r / deltaX); 
			fTheta += fDelta;
			matrix4x4 matRotZ = getMatrixRotationZ(0);
			matrix4x4 matRotX = getMatrixRotationX(0);
			matrix4x4 matRotY = getMatrixRotationY(fTheta);
			if (Math.Abs(x) >= r )
			{
				isHalfRotated = !isHalfRotated;
			}
			x += isHalfRotated ? deltaX : -deltaX;


			double z = isHalfRotated ? -Math.Sqrt(r * r - x * x) : Math.Sqrt(r * r - x * x);
			z += 4 * r;


			matrix4x4 matTrans = getMatrixTranslation(x, 1, z); ;
			

			// all the transformations that we need to do in this matrix
			matrix4x4 matWorld = createUnitMatrix();
			matWorld = matRotZ * matRotX * matRotY;
			matWorld = matWorld * matTrans;
			
			//camera
			vect3D vUp = new vect3D(0, 1, 0);
			vect3D vTarget = new vect3D(0, 0, 1);
			matrix4x4 matCameraRot = getMatrixRotationY(fYaw) * getMatrixRotationX(fXaw);
			
			vLookDir = matCameraRot * vTarget;
			vTarget = vCamera + vLookDir;
			matrix4x4 matCamera = matrixPointAt(vCamera, vTarget, vUp);
			matrix4x4 matView = matrixQuickInverse(matCamera);

			List<triangle> trianglesToRaster = new List<triangle>();

			foreach (triangle tri in meshCube.tris)
			{
				triangle triProjected = new triangle();
				triangle triTransformed = new triangle();
				triangle triViewed = new triangle();
				triTransformed.vect[0] = matWorld * tri.vect[0];
				triTransformed.vect[1] = matWorld * tri.vect[1];
				triTransformed.vect[2] = matWorld * tri.vect[2];

				vect3D normal = new vect3D();
				vect3D line1 = new vect3D();
				vect3D line2 = new vect3D();

				line1 = triTransformed.vect[1] - triTransformed.vect[0];
				line2 = triTransformed.vect[2] - triTransformed.vect[0];

				normal = line1.CrossProduct(line2);
				normal.Normalize();
			
				//projecting

				// объект виден если угол между нормалью к нему и нашей камерой меньше 90 градусов
				// будем проверять это с помощью скалярного произведения векторов
				// скалярное произведение векторов угол между которыми > 90 градусов будет < 0
				if (normal.DotProduct(triTransformed.vect[0] - vCamera) < 0.0)
				{
					// illumination (all illumination cooming from direction, not from the point)
					// свет тем ярче чем меньше угол между нормалью поверхности с направлением света
					vect3D lightDirection = new vect3D(-5.0, -5.0, -5.0);
					lightDirection.Normalize();

					double dp = normal.DotProduct(lightDirection);
					triProjected.luminosity = dp;
					triTransformed.luminosity = dp;

					triViewed.vect[0] = matView * triTransformed.vect[0];
					triViewed.vect[1] = matView * triTransformed.vect[1];
					triViewed.vect[2] = matView * triTransformed.vect[2];

					triProjected.vect[0] = matProj * triViewed.vect[0];
					triProjected.vect[1] = matProj * triViewed.vect[1];
					triProjected.vect[2] = matProj * triViewed.vect[2];

					triProjected.vect[0] = triProjected.vect[0] / triProjected.vect[0].w;
					triProjected.vect[1] = triProjected.vect[1] / triProjected.vect[1].w;
					triProjected.vect[2] = triProjected.vect[2] / triProjected.vect[2].w;

					//offset into visible mormalized space 
					vect3D vOffsetView = new vect3D(1, 1, 0);

					triProjected.vect[0] = triProjected.vect[0] + vOffsetView;
					triProjected.vect[1] = triProjected.vect[1] + vOffsetView;
					triProjected.vect[2] = triProjected.vect[2] + vOffsetView;
				
					triProjected.vect[0].x *= canvas.Width  / 2.0;
					triProjected.vect[0].y *= canvas.Height / 2.0;
					triProjected.vect[1].x *= canvas.Width  / 2.0;
					triProjected.vect[1].y *= canvas.Height / 2.0;
					triProjected.vect[2].x *= canvas.Width  / 2.0;
					triProjected.vect[2].y *= canvas.Height / 2.0;

					trianglesToRaster.Add(triProjected);
				}
			}
			// Sort triangles from back to front 
			try {
				trianglesToRaster.Sort();
			}
			catch (Exception ex)
			{

			}
			foreach (triangle t in trianglesToRaster)
			{
				DrawTriangle(t, true, true);
			}
			}
		protected override void OnKeyDown(KeyEventArgs e)
		{
			
			switch (e.Key)
			{
				//Forward
				case Key.Up:
					vCamera = vCamera + (vLookDir * 1.001); 
					break;
				//Left
				case Key.Left:
					vCamera.x -= 1.0;
					break;
				//Right
				case Key.Right:
					vCamera.x += 1.0;
					break;
				//Back		
				case Key.Down:
					vCamera = vCamera - (vLookDir * 1.001);
					break;
				// rotate up
				case Key.W:
					fXaw += 0.05;
					break;
				// rotate down
				case Key.S:
					fXaw -= 0.05;
					break;
				// turn left
				case Key.A:
					fYaw += 0.05;
					break;
				// turn right
				case Key.D:
					fYaw -= 0.05;
					break;
			}
			base.OnKeyDown(e);
		}
	}
}
