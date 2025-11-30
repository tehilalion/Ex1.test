//package assignments.Ex1;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
    /**
     * Epsilon value for numerical computation, it serves as a "close enough" threshold.
     */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /**
     * The zero polynomial function is represented as an array with a single (0) entry.
     */
    public static final double[] ZERO = {0};

    /**
     * Computes the f(x) value of the polynomial function at x.
     *
     * @param poly - polynomial function
     * @param x
     * @return f(x) - the polynomial function value at x.
     */
    public static double f(double[] poly, double x) {
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * This function should be implemented recursively.
     *
     * @param p   - the polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    /**
     * This function computes a polynomial representation from a set of 2D points on the polynom.
     * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
     * Note: this function only works for a set of points containing up to 3 points, else returns null.
     * First we check if the points that were given are valid if not, return null.
     * Then we check the length of each point.
     * If it is over three or under two return null.
     * if lx = = 2 //the length is 2 points
     * double x1 - xx [0] y1= yy[0] x2= xx[1] y2= yy[1] // helps us organize our points
     * double m = (y1-y2)/(x1-x2) //  finds the gradient (m) of our function
     * double b = y1- m*x1 // finds the interception point of the functions
     * ans = new double {m,b}
     * if lx==3 // the length is 3 points
     * double x1 - xx [0] y1= yy[0] x2= xx[1] y2= yy[1] x3= xx[2] y3= yy[2]// helps us organize our points
     * we calculate each denominator for each point
     * if the denom is 0 return null // denominator can never be 0
     * return new points
     *
     * @param xx
     * @param yy
     * @return an array of doubles representing the coefficients of the polynom.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {

            if (lx == 2) {
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double m = (y1 - y2) / (x1 - x2);
                double b = y1 - m * x1;
                ans = new double[]{m, b};
            } else if (lx == 3) {
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];
                double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                if (denom == 0) {
                    return null;
                }
                double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2));
                double B = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3));
                double C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3);

                return new double[]{A, B, C};

            }
        }
        return ans;

    }

    /**
     * Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
     * first we check if the polys are equal
     * if not we check the max length between them giving us a new length
     * int n
     * int (i = 0; i < N; i++) // checks the f(0) first and continues from there
     * we check the results
     * double dif = Math.abs(res1 - res2);// checks the abs value
     * if (dif > EPS) sees if its bigger than eps
     *if it is return false
     * if not ans
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true if p1 represents the same polynomial function as p2.
     */
    public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;
        if (p1 == p2) {
            return ans;
        }
        int l = Math.max(p1.length, p2.length);
       int N = Math.max(1, l);
      for (int i = 0; i < N; i++) {
          double x = i;
          double res1 = f(p1, x);
          double res2 = f(p2, x);
          double dif = Math.abs(res1 - res2);
          if (dif > EPS) {
              return false;
          }
      }
        return ans;
    }


	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
     * we start by saying if the length is 0 give back "0".
     * we check what kind of symbol is infront.
     * if (i>= 1) ans += "x"; // if the index is 1 or over then "x"
     *  if (i>= 2) ans += "^" +i;// if the index is over 2 then "^"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
         for (int i = poly.length-1; i >= 0; i--) {
             double x = poly[i];
             if(x!=0) {
                 if (ans.length()!=0) {
                     if (x>0) ans += " +";
                     else ans += "";
                 }
                 ans  += x;
                 if (i>= 1) ans += "x";
                 if (i>= 2) ans += "^" +i;
             }
            }
         if (ans.length() == 0) ans = "0";
        }
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double ans = x1;

        double h = f(p1, x1) - f(p2, x1);
        double mid = (x1 + x2) / 2;
        double fm = f(p1, mid) - f(p2, mid);
        if (Math.abs(fm) < eps) {
            return mid;
        }
        if (h * fm < 0) {
            return sameValue(p1, p2, x1, mid, eps);
        } else {
            return sameValue(p1, p2, x2, mid, eps);
        }
    }


	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
     * double delta = (x2 - x1)/numberOfSegments; // finds delta - the distance between each segment
     * for (int i = 0; i < numberOfSegments; i=i+1) {
     *             double x0 = x1 + delta*i; // the starting segment
     *             double x1p = x0 + delta; // the ending segment (p- for prime)
     *       double f0  = f(p, x0); // the y of the start found using the f function we already have
     *       double f1p = f(p, x1p); // the y of the end
     *       find the changing point between the two polys dx and dy
     * by using pitagoras  (Math.sqrt(dx*dx + dy*dy);) we find the distance between the points
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length (double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
       double delta = (x2 - x1)/numberOfSegments;
        for (int i = 0; i < numberOfSegments; i=i+1) {
            double x0 = x1 + delta*i;
            double x1p = x0 + delta;

           double f0  = f(p, x0);
           double f1p = f(p, x1p);

           double dx= x1p - x0;
           double dy= f1p-f0;

           ans+= Math.sqrt(dx*dx + dy*dy);
        }
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
     *
     * we start by giving our helper function "intersection". This will find us our intersection point.
     * now we move on to our area function.
     *   double delta = (x2 - x1) / numberOfTrapezoid; // finds the width of each section
     *    double x0 = x1 + delta * i; // finds the left limit
     *    double x1p = x0 + delta; // finds the right limit
     * double f0 = f(p1, x0) - f(p2, x0);
     *  double f1p = f(p1, x1p) - f(p2, x1p); // finds the function that is the difference between to points
     *  if (f0 * f1p >= 0)
     *    ans += 0.5 * (Math.abs(f0) + Math.abs(f1p)) * delta; // regular equation of finding the area of a trapizoid
     *    else - if there is a intersection point then we split the function in two, the left side and the right.
     *    we check the area of each side and put them together thus giving us the sum we were looking for.
     *
     *
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
    */


    public static double intersection(double[] p1, double[] p2, double a, double b) {
        double xr = sameValue(p1, p2, a, b, Ex1.EPS);
        if (Double.isNaN(xr) || xr <= a || xr >= b) {
            double f0 = f(p1, a) - f(p2, a);
            double f1 = f(p1, b) - f(p2, b);
            xr = a - f0 * (b - a) / (f1 - f0);
        }
        return xr;
    }

    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        if (numberOfTrapezoid <= 0 || x1 == x2) {
            return 0;
        }
        double ans = 0.0;
        double delta = (x2 - x1) / numberOfTrapezoid;
        for (int i = 0; i < numberOfTrapezoid; i++) {
            double x0 = x1 + delta * i;
            double x1p = x0 + delta;
            double f0 = f(p1, x0) - f(p2, x0);
            double f1p = f(p1, x1p) - f(p2, x1p);
            if (f0 * f1p >= 0) {

                ans += 0.5 * (Math.abs(f0) + Math.abs(f1p)) * delta;
            }
            else {

                double xr = intersection(p1, p2, x0, x1p);
                double delta1 = xr - x0;
                double delta2 = x1p - xr;
                ans += 0.5 * Math.abs(f0) * delta1;
                ans += 0.5 * Math.abs(f1p) * delta2;
            }
        }
        return ans;
    }



	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * we start with getA function gets a monom and changes it to be a coefficient number.
     * if (!monom.contains("x")) {
     *  return Double.parseDouble(monom);// if it doesnt contain x skip it.
     *  but if it does -int xi = monom.indexOf("x"); // finds the index x is on
     *  if (monom.substring(0, xi).equals("-")) {
     *    return -1.0;// if the number starts  with - return -1.
     * We continue to function getB which helps us find the power of the function.
     * if the function doesnt have x skip it if it contains "^" find the index its on and add 1
     * (all of these functions connect to our string from poly function)
     * Now we move on to out getPloyFromString function
     * if (p == null || p.length() == 0) return ZERO; // return 0 if the length is 0 or null
     *     p = p.replaceAll(" ", "");
     *         p = p.replace("-", "+-");
     *         if (p.startsWith("+")) p = p.substring(1); // replaces our symbols
     *  String[] terms = p.split("\\+"); // splits pur function so we can calculate
     *  maxpower finds the highest power there is.
     *  once we do this we add 1 to find the correct index.
     *  for (String t : terms) {
     *             if (t.length() == 0) continue; // we move on to our monoms and checks them according to our functions.
     *  double a = getA(t);
     *             int b = getB(t);
     *             ans[b] += a;
     *             return ans;
     *
	 * @param p - a String representing polynomial function.
	 * @return
	 */
    public static double getA(String monom) {
        monom = monom.replaceAll(" ", "");
        if (!monom.contains("x")) {
            return Double.parseDouble(monom);
        }
        int xi = monom.indexOf("x");
        if (xi == 0) {
            return 1.0;
        }
        if (monom.substring(0, xi).equals("-")) {
            return -1.0;
        }
        return Double.parseDouble(monom.substring(0, xi));
    }

    public static int getB(String monom) {
        monom = monom.replaceAll(" ", "");
        if (!monom.contains("x")) return 0;
        if (monom.contains("^")) {
            return Integer.parseInt(monom.substring(monom.indexOf("^") + 1));
        }
        return 1;
    }


    public static double[] getPolynomFromString(String p) {
        if (p == null || p.length() == 0) return ZERO;
        p = p.replaceAll(" ", "");
        p = p.replace("-", "+-");
        if (p.startsWith("+")) p = p.substring(1);

        String[] terms = p.split("\\+");

        int maxPower = 0;
        for (String t : terms) {
            if (t.length() == 0) continue;
            int b = getB(t);
            if (b > maxPower) maxPower = b;
        }
        double[] ans = new double[maxPower + 1];

        for (String t : terms) {
            if (t.length() == 0) continue;

            double a = getA(t);
            int b = getB(t);
            ans[b] += a;
        }
        return ans;
    }


	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     * we first check the max length of the polys and make a new length according to the max result we got.
     * we start a loop
     * double a = 0
     * double b = 0
     * if i<pi.length a=p1[i] // if the index is smaller than the length a = the index
     * same for b and p2
     * when a and b have the same index add them together.
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
        int l = Math.max (p1.length, p2.length);
      double [] ans = new double [l];
        for (int i = 0; i < l; i++) {
            double a = 0;
            double b = 0;
            if (i < p1.length) {
                a = p1[i];
            }
            if (i < p2.length) {
                b = p2[i];
            }
            ans[i] = a + b;
        }
        return ans;
    }

	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     * we add the lengthes and minus 1 in order to get the right index this is our new length.
     * we start with sum being 0
     *
	 * @param p1
	 * @param p2
	 * @return
     *
	 */
	public static double[] mul(double[] p1, double[] p2) {
            int l= (p1.length+p2.length-1);
            double [] ans = new double [l];
            for (int i = 0; i < l; i=i+1) {
                double sum = 0.0;

             for (int  a=0; a<p1.length; a=a+1) {
                  int b = i - a;

                if (b>=0 && b<p2.length) {
                    sum += p1[a]*p2[b];
                }
             }
                ans[i]= sum;
    }
        return ans;
        }
	/**
	 * This function computes the derivative of the p0 polynomial function.
     * The method returns a new array representing the derivative polynomial.
     *  * If the input polynomial is  null or has length 1 (a constant polynomial), the method returns zero
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
     if(po!=null && po.length>1) {
     int len = po.length;
     ans = new double [len-1];
     for (int i = 0; i < ans.length; i=i+1) {
     ans[i] = po[i+1] * (i+1);
        }
     }
     return ans;
    }





}
