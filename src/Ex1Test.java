//package assignments.Ex1;
import com.sun.tools.javac.Main;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 * @author boaz.ben-moshe
 */

class Ex1Test {
    static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1,0.1,3};
    static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};
    static double[] po3 = {2,1,-0.7, -0.02,0.02};
    static double[] po4 = {-3, 0.61, 0.2};

    @Test
    /**
     * this test checks if the function will return null
     */
    public void testPolynomFromPoints() {
        double[] xx = {1.0, 3.0};
        double[] yy = {3.0, 7.0};
        double[] result = Ex1.PolynomFromPoints(xx, yy);
        assertNotNull(result, "result should not be null");

        /**
         * this function checks the length
         */
        assertEquals(2, result.length, "result length should be 2");

        /**
         * this function checks if the gradient and interception point are correct
         */
        double expectedm = 2.0;
        double expectedb = 1.0;
        assertEquals(expectedm, result[0], 0.001, "gradient should be correct");
        assertEquals(expectedb, result[1], 0.001, "intercept should be correct");
    }






    @Test
    /**
     * Tests that f(x) == poly(x).
     */
    void testF() {
        double fx0 = Ex1.f(po1, 0);
        double fx1 = Ex1.f(po1, 1);
        double fx2 = Ex1.f(po1, 2);
        assertEquals(fx0, 2, Ex1.EPS);
        assertEquals(fx1, 4, Ex1.EPS);
        assertEquals(fx2, 6, Ex1.EPS);
    }
    @Test
    /**
     * Tests that p1(x) + p2(x) == (p1+p2)(x)
     */
    void testF2() {
        double x = Math.PI;
        double[] po12 = Ex1.add(po1, po2);
        double f1x = Ex1.f(po1, x);
        double f2x = Ex1.f(po2, x);
        double f12x = Ex1.f(po12, x);
        assertEquals(f1x + f2x, f12x, Ex1.EPS);
    }
    @Test
    /**
     * Tests that p1+p2+ (-1*p2) == p1
     */
    void testAdd() {
        double[] p12 = Ex1.add(po1, po2);
        double[] minus1 = {-1};
        double[] pp2 = Ex1.mul(po2, minus1);
        double[] p1 = Ex1.add(p12, pp2);
        assertTrue(Ex1.equals(p1, po1));
    }
    @Test
    /**
     * Tests that p1+p2 == p2+p1
     */
    void testAdd2() {
        double[] p12 = Ex1.add(po1, po2);
        double[] p21 = Ex1.add(po2, po1);
        assertTrue(Ex1.equals(p12, p21));
    }
    @Test
    /**
     * Tests that p1+0 == p1
     */
    void testAdd3() {
        double[] p1 = Ex1.add(po1, Ex1.ZERO);
        assertTrue(Ex1.equals(p1, po1));
    }
    @Test
    /**
     * Tests that p1*0 == 0
     */
    void testMul1() {
        double[] p1 = Ex1.mul(po1, Ex1.ZERO);
        assertTrue(Ex1.equals(p1, Ex1.ZERO));
    }
    @Test
    /**
     * Tests that p1*p2 == p2*p1
     */
    void testMul2() {
        double[] p12 = Ex1.mul(po1, po2);
        double[] p21 = Ex1.mul(po2, po1);
        assertTrue(Ex1.equals(p12, p21));
    }
    @Test
    /**
     * Tests that p1(x) * p2(x) = (p1*p2)(x),
     */
    void testMulDoubleArrayDoubleArray() {
        double[] xx = {0,1,2,3,4.1,-15.2222};
        double[] p12 = Ex1.mul(po1, po2);
        for(int i = 0;i<xx.length;i=i+1) {
            double x = xx[i];
            double f1x = Ex1.f(po1, x);
            double f2x = Ex1.f(po2, x);
            double f12x = Ex1.f(p12, x);
            assertEquals(f12x, f1x*f2x, Ex1.EPS);
        }
    }
    @Test
    /**
     * Tests a simple derivative examples - till ZERO.
     */
    void testDerivativeArrayDoubleArray() {
        double[] p = {1,2,3}; // 3X^2+2x+1
        double[] pt = {2,6}; // 6x+2
        double[] dp1 = Ex1.derivative(p); // 2x + 6
        double[] dp2 = Ex1.derivative(dp1); // 2
        double[] dp3 = Ex1.derivative(dp2); // 0
        double[] dp4 = Ex1.derivative(dp3); // 0
        assertTrue(Ex1.equals(dp1, pt));
        assertTrue(Ex1.equals(Ex1.ZERO, dp3));
        assertTrue(Ex1.equals(dp4, dp3));
    }
    @Test
    /**
     * Tests the parsing of a polynom in a String like form.
     */
    public void testFromString() {
        double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
        String sp2 = "3.1x^2 +2.3x -1.1";
        String sp = Ex1.poly(p);
        double[] p1 = Ex1.getPolynomFromString(sp);
        double[] p2 = Ex1.getPolynomFromString(sp2);
        boolean isSame1 = Ex1.equals(p1, p);
        boolean isSame2 = Ex1.equals(p2, p);
        if(!isSame1) {fail();}
        if(!isSame2) {fail();}
        assertEquals(sp, Ex1.poly(p1));
    }
    @Test
    /**
     * Tests the equality of pairs of arrays.
     */
    public void testEquals() {
        double[][] d1 = {{0}, {1}, {1,2,0,0}};
        double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
        double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
        for(int i=0;i<d1.length;i=i+1) {
            assertTrue(Ex1.equals(d1[i], d2[i]));
        }
        for(int i=0;i<d1.length;i=i+1) {
            assertFalse(Ex1.equals(d1[i], xx[i]));
        }
    }

    @Test
    /**
     * Tests is the sameValue function is symmetric.
     */
    public void testSameValue2() {
        double x1=-4, x2=0;
        double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
        double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
        assertEquals(rs1,rs2, Ex1.EPS);
    }
    @Test
    /**
     * Test the area function - it should be symmetric.
     */
    public void testArea() {
        double x1=-4, x2=0;
        double a1 = Ex1.area(po1, po2, x1, x2, 100);
        double a2 = Ex1.area(po2, po1, x1, x2, 100);
        assertEquals(a1,a2, Ex1.EPS);
    }
    @Test
    /**
     * Test the area f1(x)=0, f2(x)=x;
     */
    public void testArea2() {
        double[] po_a = Ex1.ZERO;
        double[] po_b = {0,1};
        double x1 = -1;
        double x2 = 2;
        double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
        double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
        double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
        double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
        double area =2.5;
        assertEquals(a1,area, Ex1.EPS);
        assertEquals(a2,area, Ex1.EPS);
        assertEquals(a3,area, Ex1.EPS);
        assertEquals(a100,area, Ex1.EPS);
    }
    @Test
    /**
     * Test the area function.
     */
    public void testArea3() {
        double[] po_a = {2,1,-0.7, -0.02,0.02};
        double[] po_b = {6, 0.1, -0.2};
        double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
        double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
        double area = 58.5658;
        assertEquals(a1,area, Ex1.EPS);
    }

 /** @Test
   public void testGetB1() {}
    String [] monoms = {"x" , "1x^1", "3", "-32.1", "-3.1x^3"};
    int [] res= {1,1,0,0,3};
    for(int i=0; i<res.length; i=i+1) {
        int r= Ex1.getB(monoms[i]);
        assertEquals(res[i], r);
    }
*/
        @Test
        /**
         * checks when its a normal poly
         */
        public void testNormalPolynomial() {
            double[] poly = {2, 0, 3.1, -1.2};
            String expected = "-1.2x^3 +3.1x^2 +2.0";
            assertEquals(expected, Ex1.poly(poly));
        }

        @Test
        /**
         * checks when a poly is all Zero
         */
        public void testAllZeros() {
            double[] poly = {0, 0, 0};
            String expected = "0";
            assertEquals(expected, Ex1.poly(poly));
        }

        @Test
        /**
         * checks when a poly is with only one term
         */
        public void testSingleTerm() {
            double[] poly = {0, 0, 5};
            String expected = "5.0x^2";
            assertEquals(expected, Ex1.poly(poly));
        }

        @Test
        /**
         * checks when a poly has a negative
         */
        public void testNegativeCoefficients() {
            double[] poly = {-3, 0, -2};
            String expected = "-2.0x^2-3.0";
            assertEquals(expected, Ex1.poly(poly));
        }

        @Test
        /**
         * checks when poly is empty
         */
        public void testEmptyArray() {
            double[] poly = {};
            String expected = "0";
            assertEquals(expected, Ex1.poly(poly));
        }

    @Test
    /**
     * checks what happens when its a straight line
     */
    public void testLength1() {
        double[] p = {0, 1};
        double x1 = 0;
        double x2 = 3;
        int segments = 1000;
        double result = Ex1.length(p, x1, x2, segments);
        double expected = (x2 - x1) * Math.sqrt(2);
        assertEquals(expected, result - x1, 0.01,
                "Length of straight line y=x should be √2 * (x2−x1)");
    }

@Test
/**
 * checks when it's a horizontial line
 */
    public void testLength2() {
        double[] p = {5};
        double x1 = -2;
        double x2 = 4;
        int segments = 1000;
        double result = Ex1.length(p, x1, x2, segments);
        double expected = (x2 - x1);
        assertEquals (expected, result-x1, 0.01);
}
    @Test
    /**
     * checks when poly is zero
     */
    public void testLength3() {
        double[] p = {0}; // f(x) = 0
        double x1 = 1;
        double x2 = 5;
        int segments = 1000;
        double result = Ex1.length(p, x1, x2, segments);
        double expected = (x2 - x1);
        assertEquals(expected, result - x1, 0.01);
    }
    @Test
    /**
     * checks when its a parabola
     */
    public void testLength4() {
        double[] p = {0, 0, 1}; // f(x) = x^2
        double x1 = 0;
        double x2 = 1;
        int segments = 2000;
        double result = Ex1.length(p, x1, x2, segments);
        double expected = 0.25 * (2 * Math.sqrt(5) + Math.log(2 + Math.sqrt(5)));
        assertEquals(expected, result - x1, 0.01);
    }









}
