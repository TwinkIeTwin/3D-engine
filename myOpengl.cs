using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _3D_engine
{
	class MyOpengl
	{

		// вектор показывает относительно направление к заданной точки от стартовой (0, 0)
		class Vector
		{
			private double[] Coordinates { get; set; }

			public double Length { get; set; }

			public int Size {get; set;}

			public Vector(int size)
			{
				this.Size = size;
				this.Coordinates = new double[size];
			}
			// будем пересчитывать длинну вектора после изменения в его координатах
			private void ReCalcLength()
			{
				double len = 0;
				foreach (var i in Coordinates)
				{
					len += i * i;
				}
				len = Math.Sqrt(len);
			}

			public Vector(int size, params double[] coordinates)
			{
				this.Size = size;
				this.Coordinates = coordinates;
				ReCalcLength();
			}

			// будем получать направление с углом по середине между двумя заданными
			public static Vector operator +(Vector v1, Vector v2)
			{
				Vector vResult = new Vector(v1.Size);
				int coordIndex = 0;
				foreach (var i in v1.Coordinates)
				{
					foreach (var j in v2.Coordinates)
					{
						vResult.Coordinates[coordIndex++] = i + j;
					}
				}
				vResult.ReCalcLength();
				return vResult;
			}

			
			// чтобы получить направление из v1 в v2: v2 - v1
			public static Vector operator -(Vector v1, Vector v2)
			{
				Vector vResult = new Vector(v1.Size);
				int coordIndex = 0;
				foreach (var i in v1.Coordinates)
				{
					foreach (var j in v2.Coordinates)
					{
						vResult.Coordinates[coordIndex++] = i - j;
					}
				}
				vResult.ReCalcLength();
				return vResult;
			}

			public static Vector operator +(Vector v, double c)
			{
				Vector vResult = new Vector(v.Size);
				int vSize = v.Coordinates.Length;
				for (int i = 0; i < vSize; i++)
				{
					vResult.Coordinates[i] = v.Coordinates[i] + c;
				}
				vResult.ReCalcLength();
				return vResult;
			}


			public static Vector operator -(Vector v, double c)
			{
				Vector vResult = new Vector(v.Size);
				int vSize = v.Coordinates.Length;
				for (int i = 0; i < vSize; i++)
				{
					vResult.Coordinates[i] = v.Coordinates[i] - c;
				}
				vResult.ReCalcLength();
				return vResult;
			}
			// удлинить вектор в c раз не изменив направление
			public static Vector operator *(Vector v, double c)
			{
				Vector vResult = new Vector(v.Size);
				int vSize = v.Coordinates.Length;
				for (int i = 0; i < vSize; i++)
				{
					vResult.Coordinates[i] = v.Coordinates[i] * c;
				}
				vResult.ReCalcLength();
				return vResult;
			}

			public static Vector operator /(Vector v, double c)
			{
				Vector vResult = new Vector(v.Size);
				int vSize = v.Coordinates.Length;
				for (int i = 0; i < vSize; i++)
				{
					vResult.Coordinates[i] = v.Coordinates[i] / c;
				}
				vResult.ReCalcLength();
				return vResult;
			}

			public Vector GetDirectionTo(Vector vTarget)
			{
				return vTarget - this;
			}

			public Vector GetDirectionFromTo(Vector vFrom, Vector vTo)
			{
				return vTo - vFrom;
			}

			// сделать вектор единичной длинны
			public void Normalize()
			{
				if (this.Length == 0)
				{
					throw new Exception("Длинна вектора была 0 при попытке нормализации");
				}
				for (int i = 0; i < this.Coordinates.Length; i++)
				{
					this.Coordinates[i] /= this.Length;
				}
			}

		}		
		
		class Vector2D
		{
			public double X { get; set;}
			public double Y { get; set;}

			public Vector2D(double x, double y)
			{
				this.X = x;
				this.Y = y;
			}
		

			public static Vector2D operator +(Vector2D v1, Vector2D v2)
			{
				return new Vector2D(v1.X + v2.X, v1.Y + v2.Y);
			}

		}

		class Vector3D
		{
			

		}


	}
}
