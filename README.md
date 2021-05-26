# A simple Analytic Element Method (AEM) toolbox

This repository contains the two Python3 toolboxes and a simple tutorial for the accompanying paper in [Water Resources Research](https://agupubs.onlinelibrary.wiley.com/journal/19447973) (current status: under review).

## Installation

To use these toolboxes, simply download the Python files [`toolbox_AEM.py`](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_AEM.py) and [`toolbox_MCMC.py`](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/toolbox_MCMC.py) (if desired) and copy them into your working directory. The MCMC toolbox is optional and only required if you also wish to use my MCMC implementation, the AEM toolbox can work as a standalone. Both toolboxes can be imported into your Python script by including the following code snippet:

```
from toolbox_AEM import *
from toolbox_MCMC import *
```

And that's it! 

## Tutorials


I have provided a few simple tutorials to help you familiarize yourself with the toolbox. The [**basic tutorial**](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/Tutorial%2001%20Basic%20AEM) covers the fundamentals of constructing a deterministic flow model with this toolbox, and is available as both a [Python file](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2001%20Basic%20AEM/basic_tutorial.py) and a [Jupyter Notebook](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2001%20Basic%20AEM/basic_tutorial.ipynb). 

The [**advanced tutorial**](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/Tutorial%2002%20Uncertainty%20Estimation) expands the scope of the basic tutorial by showing how to prepare the Analytic Element Model for uncertainty quantification. This tutorial is also available as both a [Python file](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2002%20Uncertainty%20Estimation/uncertainty_estimation_example.py) and a [Jupyter Notebook](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/blob/main/Tutorials/Tutorial%2002%20Uncertainty%20Estimation/uncertainty_estimation_example.ipynb).

For users with an interest in reproducing the some or all of the results in accompanying manuscript, I have also [uploaded the Python files](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Manuscript%20files) required to reproduce the figures in the main manuscript and the supporting information.



<!DOCTYPE html>
<html>

  <head>
    <meta name="description" content="math.js | basic usage">
    <title>math.js | basic usage</title>
    <script src="https://unpkg.com/mathjs/lib/browser/math.js"></script>
    
    <script src="https://d3js.org/d3.v4.min.js"></script>
    <script src="https://d3js.org/d3-hsv.v0.1.min.js"></script>
    <script src="https://d3js.org/d3-contour.v1.min.js"></script>
    <script src="https://d3js.org/d3-scale-chromatic.v1.min.js"></script>

    <style>

    </style>
  </head>

<!-- Create a div where the graph will take place -->
<div id="my_dataviz">
  <svg id="click" xmlns="http://www.w3.org/2000/svg">
      <defs>
          <g id="pointer" transform="scale(0.5)">
              <circle cx="0" cy="0" r="20" id="dragcircle" />
          </g>
      </defs>
  </svg>
</div>



  <body style='overflow:hidden'>
  
    
    <script>
    
      const vw = Math.max(document.documentElement.clientWidth || 0, window.innerWidth || 0)
      const vh = Math.max(document.documentElement.clientHeight || 0, window.innerHeight || 0)
      
      var height 	= math.min([vw,vh]); // 450
      var width 	= math.min([vw,vh]);
      
      var height_offset = 0.05*height
      var width_offset = 0.05*width
    
      var p1 				= [-0.5, -0.5];
      var p2 				= [0.5, 0.5];
      
      var w1 				= [0.2, -0.2];
      var rw1 			= 0.01*height/450; // Well screen radius
      var sw1 			= -1; // Well strength
      
      var w2 				= [-0.2, 0.2];
      var rw2 			= 0.01*height/450; // Well screen radius
      var sw2 			= 1; // Well strength
      
      var vres 			= 7;
      
      var dragging  = false;
      
      var l1 			 	= math.complex(
      	math.cos(math.dotMultiply(-0.25,math.PI)),
        math.sin(math.dotMultiply(-0.25,math.PI)))
      var l2 			 	= math.complex(
      	math.cos(math.dotMultiply(0.25,math.PI)),
        math.sin(math.dotMultiply(0.25,math.PI)))
      var l3 			 	= math.complex(
      	math.cos(math.dotMultiply(0.75,math.PI)),
        math.sin(math.dotMultiply(0.75,math.PI)))

      var z1 = angle_to_complex([math.dotMultiply(-0.4,math.PI)])[0];
      var z2 = angle_to_complex([math.dotMultiply(0.1,math.PI)])[0];
      var z3 = angle_to_complex([math.dotMultiply(0.4,math.PI)])[0];
      
      
      [a,b,c,d] = get_moebius_coefficients(z1,z2,z3,l1,l2,l3)
      
			var z4 = math.complex(
      	math.cos(math.dotMultiply(-0.75,math.PI)),
        math.sin(math.dotMultiply(-0.75,math.PI)))
        
      var p4 = moebius(z4,a,b,c,d);
      
      // print([a,b,d,c])
      
      var pa 	= a;
      var pb 	= b;
      var pc 	= c;
      var pd 	= d;
      
      var vertices 	= linspace_complex(p1, p2, vres);
      var zc 			 	= linspace_complex_segments(p1, p2, vres);
      var nvec 			= normal_vector(p1,p2);
      var A 				= noflowboundary_gradients(vertices,zc,nvec);
      var solvec 		= solution_vector(zc,nvec);
      var s 				= solve_linear_system(A,solvec);
      var x_res			= 41;
      var y_res 		= 41;
      var grid 			= construct_grid([-1,1],[-1,1],x_res,y_res);
      var output 		= evaluate(grid, vertices, s);


      // helper function to output formatted results.
      function print(value) {
        var precision = 14;
        document.write(math.format(value, precision) + '<br>');
      }
      
      function angle_to_complex(angle) {
      	var results 	= [];
        for (i = 0; i < angle.length; i ++) {
        	results.push(math.complex(
            math.cos(angle[i]),
            math.sin(angle[i]) ) )
        }
        return results;
      }

      // This function creates 11 points equally spaced between a tuple p1 and a tuple p2
      function linspace(start,end,resolution) {
        var spacing = [];
        // Go through a for-loop
        var i;
        for (i = 0; i < resolution; i++) {
        	spacing.push(start + (end-start)*i/(resolution-1))
        }
        return spacing; // The function returns the linspace between p1 and p2
      }

      // This function creates 11 points equally spaced between a tuple p1 and a tuple p2
      function linspace_complex(p1, p2, resolution) {
        var spacing = [];
        // Go through a for-loop
        var i;
        for (i = 0; i < resolution; i++) {
        	spacing.push(math.complex(
          	p1[0] + (p2[0] -p1[0] )*i/(resolution-1),
            p1[1]  + (p2[1]-p1[1])*i/(resolution-1)))
        }
        return spacing; // The function returns the linspace between p1 and p2
      }

      // This function creates 11 points equally spaced between a tuple p1 and a tuple p2
      function linspace_complex_segments(p1, p2, resolution) {
        var spacing = [];
        // Go through a for-loop
        xdif = (p2[0]-p1[0])/(resolution-1)/2
        ydif = (p2[1]-p1[1])/(resolution-1)/2
        var i;
        for (i = 0; i < resolution-1; i++) {
        	spacing.push(math.complex(
          	p1[0] + (p2[0] -p1[0] )*i/(resolution-1) + xdif,
            p1[1]  + (p2[1]-p1[1])*i/(resolution-1) + ydif))
        }
        return spacing; // The function returns the linspace between p1 and p2
      }
      
			// =============================================================================
      // EVALUATE THE MODEL
      // =============================================================================
      
      function evaluate(x, vertices, s) {
      
      	// Segments are between the vertices, so we have one less than vertices
        var segments = vertices.length-1
        
        // Go through a for-loop
        var i;
        for (i = 0; i < segments; i++) {
        
        	// Project the control points on the real line between -1 and 1
          Z = cdivide(
          	math.subtract(
            	math.dotMultiply(x, 2),
              cadd(vertices[i],vertices[i+1])),
            csub(vertices[i+1],vertices[i]));
              
          // Calculate term 1
          term1 	= math.dotMultiply(
          	math.add(
            	Z,
              1),
           	math.log(
            	math.dotDivide(
                math.subtract(
                  Z,
                  1),
                math.add(
                  Z,
                  1)) ) )
          
          // Calculate term 2
          term2 	= math.dotMultiply(
          	math.subtract(
            	Z,
              1),
           	math.log(
            	math.dotDivide(
                math.subtract(
                  Z,
                  1),
                math.add(
                  Z,
                  1))))
          
          numer = math.dotMultiply(
              s[i][0],
              math.subtract(
                  term1,
                  term2))
          
          denom = math.dotMultiply(
          	4*math.PI,
            math.complex(0,1) )
          
          temp 	= math.dotDivide(
          	numer,
            denom)

          // Add the results to the output vector
          if (i == 0) {
          	// In the first iteration, initialize the output vector
          	var results = temp
          } else {
          	results 		 = math.add(
            	results,
              temp)
          }
        }
        
        // Add the effect of the well
        
        //# Find the indices of wells      
        //dist            = np.abs(z-self.zc)
        //idx_inside      = np.where(dist < self.rw)[0]
        
        //# Correct the coordinates ---------------------------------------------
        //# Center the evaluation points on the well
        //zs  = z.copy()-self.zc
        
        //# Snap points inside the well to the well edge
        //zs[idx_inside] = self.rw + 0j    
        
        //# Calculate the complex potential
        //res = -self.strength/(2*np.pi)*np.log(zs) - temp
        
				
        // Get the local coordinates
        zs 	= math.subtract(
          x,
          math.complex(w1[0],w1[1]))
          
        //print(zs)

        // Truncate any coordinates inside the well screen
        for (j = 0; j < zs.length; j++) {
          radius = math.abs(zs[j])
          if (radius < rw1) {
            zs[j] = math.dotMultiply(zs[j],rw1/radius);
          }}

        // Calculate the gradient denominator
        temp = math.dotDivide(
        	math.dotMultiply(
          	sw1,
            math.log(zs)),
         	math.dotMultiply(
          	2,
            math.PI))
        

        // Add this to the stack (actually subtract a negative value)
        results = math.subtract(
          results,
          temp)
          
          
          
        // Get the local coordinates
        zs 	= math.subtract(
          x,
          math.complex(w2[0],w2[1]))
          
        //print(zs)

        // Truncate any coordinates inside the well screen
        for (j = 0; j < zs.length; j++) {
          radius = math.abs(zs[j])
          if (radius < rw2) {
            zs[j] = math.dotMultiply(zs[j],rw2/radius);
          }}

        // Calculate the gradient denominator
        temp = math.dotDivide(
        	math.dotMultiply(
          	sw2,
            math.log(zs)),
         	math.dotMultiply(
          	2,
            math.PI))
        

        // Add this to the stack (actually subtract a negative value)
        results = math.subtract(
          results,
          temp)
          
          
          
          
				//print(x)

				// Also add the Möbius flow
        results 			= math.add(
        	results,
          moebius_inverse(
            x,
            a,
            b,
            c,
            d))
            

           
			
				
        //// Also add the uniform flow
        //results 		 = math.add(
        //  results,
        //  x)
          
        return creal(results);
      }
     
      
			// =============================================================================
      // CONSTRUCT GRID
      // =============================================================================
      
      function construct_grid(x_limits,y_limits,x_res,y_res) {
      
      	// This function creates a rectangular grid
        var grid 	= [];
        var ydif  = y_limits[1]-y_limits[0];
        var xdif 	= x_limits[1]-x_limits[0];
        
        // Go through a for-loop
        var i;
        var j;
        for (j = 0; j < y_res; j++) {
          for (i = 0; i < x_res; i++) {
						grid.push(math.complex(
            	x_limits[0] + xdif*i/(x_res-1)+0.0000000001,
              y_limits[0] + ydif*j/(y_res-1)+0.0000000001));
          }
        }
        
				//print(grid)
        return grid; 
      }
      
      
			// =============================================================================
      // SOLVE SYSTEM OF EQUATIONS
      // =============================================================================
      
      function solve_linear_system(A,c) {
      
      	// In this case, the solution vector is really simple: Just standard uniform
        // flow from east to west with a gradient of one. As such:

        return math.lusolve(A, c); 
      }
      
      
			// =============================================================================
      // FIND NORMAL VECTOR
      // =============================================================================
      
      function normal_vector(p1,p2) {
      
      	// Define the vector between the endpoints, get its orthogonal part
      	nvec 	= [
        	p1[1] - p2[1],
          p2[0] - p1[0]];

        // Normalize this vector
        nvec = math.dotDivide(
        	nvec,
          math.sqrt(math.pow(nvec[0],2) + math.pow(nvec[1],2)));

        return nvec;
      }
      
			// =============================================================================
      // SOLUTION VECTOR
      // =============================================================================
      
      function solution_vector(zc, nvec) {
      
      	// In this case, the solution vector is really simple: Just standard uniform
        // flow from east to west with a gradient of one. As such:
        gradient 	= 1;
        solvec 		= [];
        
        // Go through a for-loop
        var i;
        for (i = 0; i < zc.length; i++) {
        

          
          // Initialize the background gradient
          //temp 	= math.complex(-gradient,0);
          //val 	= moebius_inverse(zc[i],a,b,c,d)
          
          //temp 	= cdivide(
          //	math.subtract(
          //		math.dotMultiply(a,d),
          //    math.dotMultiply(b,c)),
          //  math.pow(
          //  	math.subtract(
          //    	math.dotMultiply(
          //      	val,
          //        c),
          //      a),
          //    2))
          
          temp 	= cdivide(
          	math.subtract(
          		math.dotMultiply(a,d),
              math.dotMultiply(b,c)),
            math.pow(
            	math.subtract(
              	a,
              	math.dotMultiply(
                	zc[i],
                  c)),
              2))
          
              
          // Convert it to the partial derivatives of the hydraulic potential
          temp = math.conj(temp)
          
          // We need to compensate the OPPOSITE of the induced gradient
          temp = math.dotMultiply(temp,math.complex(-1,0))
          
          
          // Add the well gradient
          // grad[idx_valid] = -self.strength/(2*np.pi)/zs[idx_valid]
          
          // Get the local coordinates
          zs 	= math.subtract(
            zc[i],
            math.complex(w1[0],w1[1]))
            
          // Truncate any coordinates inside the well screen
          radius = math.abs(zs)
          if (radius < rw1) {
            zs 	= math.dotMultiply(
            	zs,rw1/radius);
          }
          
          // Calculate the gradient denominator
          denom = math.dotMultiply(
          	2,
            math.dotMultiply(
            	math.PI,
              zs))
          
          // Calculate the gradient numerator
          numer = sw1
          
          // Calculate the well's induced complex potential gradient
          well_effect = math.dotDivide(
            	numer,
              denom)
              
          // Convert it to the partial derivatives of the hydraulic potential
          well_effect = math.conj(well_effect)
              
          // Add this to the stack (actually subtract a negative value)
          temp 	= math.add(
          	temp,
            well_effect)
            
            
            
            
          // Get the local coordinates
          zs 	= math.subtract(
            zc[i],
            math.complex(w2[0],w2[1]))
            
          // Truncate any coordinates inside the well screen
          radius = math.abs(zs)
          if (radius < rw2) {
            zs 	= math.dotMultiply(
            	zs,rw2/radius);
          }
          
          // Calculate the gradient denominator
          denom = math.dotMultiply(
          	2,
            math.dotMultiply(
            	math.PI,
              zs))
          
          // Calculate the gradient numerator
          numer = sw2
          
          // Calculate the well's induced complex potential gradient
          well_effect = math.dotDivide(
            	numer,
              denom)
              
          // Convert it to the partial derivatives of the hydraulic potential
          well_effect = math.conj(well_effect)
              
          // Add this to the stack (actually subtract a negative value)
          temp 	= math.add(
          	temp,
            well_effect)
            
            
            
            
            
            

          
          // Now calculate the inner product between the partial derivatives 
          // and the normal vector
          

          temp = math.add(
            math.dotMultiply(
              nvec[0],
              temp.re),
            math.dotMultiply(
              nvec[1],
              temp.im));
              
          
          
          // Append the results to the matrix
					solvec.push(temp);
          
        } 
        
        return solvec; // The function returns the linspace between p1 and p2
      }

			// =============================================================================
      // INFLUENCE MATRIX
      // =============================================================================
      
      function noflowboundary_gradients(vertices,zc, nvec) {
      
      	// Segments are between the vertices, so we have one less than vertices
        var segments = vertices.length-1

				// Initialize the output array
        var results = [];
        
        // Go through a for-loop
        var i;
        for (i = 0; i < segments; i++) {
        
        	// Project the control points on the real line between -1 and 1
          var Z = cdivide(
          	math.subtract(
            	math.dotMultiply(zc, 2),
              cadd(vertices[i],vertices[i+1])),
            csub(vertices[i+1],vertices[i]));

          // Convert to dOmega(Z)/dZ
          temp = math.dotDivide(
            math.complex(0,1),
            math.subtract(
            	math.PI,
              math.dotMultiply(
              	math.PI,
                cpow(
                	Z,
                  2))));
                  
          // Multiply with dZ/dz to obtain dOmega(Z)/dz
          temp = math.dotDivide(
            	math.dotMultiply(
          			temp,
                2),
              csub(vertices[i+1],vertices[i]) )
          
					// We require the partial derivatives of the real component, get the conjugate
          temp = math.conj(temp)
          
          // Now calculate the inner product between the partial derivatives 
          // and the normal vector
          temp = math.add(
          	math.dotMultiply(
            	nvec[0],
              creal(temp)),
            math.dotMultiply(
            	nvec[1],
              cimag(temp)));
              
          // Append the results to the matrix
          results.push(temp);
        } 

        return results; // Take its transpose: math.transpose(results)
      }
      
			// =============================================================================
      // FUNCTION MOEBIUS
      // =============================================================================
      
      function moebius(z, a, b, c, d) {
        
        results 	= cdivide(
          math.add(
            math.dotMultiply(z,a),
            b),
          math.add(
            math.dotMultiply(z,c),
            d) )
            
         return results;
      }
      
      function moebius_inverse(z, a, b, c, d) {
        
        results 	= math.dotDivide(
          math.subtract(
            b,
            math.dotMultiply(z,d)),
          math.subtract(
            math.dotMultiply(z,c),
            a) )
            
         return results;
      }
      
      function get_moebius_coefficients(z1,z2,z3,l1,l2,l3) {
				
        // var l1 = math.complex(1,1);
        // var l2 = math.complex(1,1);
        // var l3 = math.complex(1,1);
        
        var a = math.det(
        	[[cmultiply(l1,z1), z1, 1],
           [cmultiply(l2,z2), z2, 1],
           [cmultiply(l3,z3), z3, 1]])
        
        var b = math.det(
        	[[cmultiply(l1,z1), l1, z1],
           [cmultiply(l2,z2), l2, z2],
           [cmultiply(l3,z3), l3, z3]])
        
        var c = math.det(
        	[[l1, z1, 1],
           [l2, z2, 1],
           [l3, z3, 1]])
        
        var d = math.det(
        	[[cmultiply(l1,z1), l1, 1],
           [cmultiply(l2,z2), l2, 1],
           [cmultiply(l3,z3), l3, 1]])
        

        

        
        //mat1 	= [
        //	[math.subtract(z2,z3), 	
        //   math.subtract(
        //   	math.dotMultiply(z1,z3),
        //    math.dotMultiply(z1,z2))],
        //  [math.subtract(z2,z1), 	
        //   math.subtract(
        //   	math.dotMultiply(z1,z3),
        //    math.dotMultiply(z2,z3))]];
      
        //mat2 	= [
        //	[math.subtract(l2,l3), 	
        //   math.subtract(
        //   	math.dotMultiply(l1,l3),
        //    math.dotMultiply(l1,l2))],
        //  [math.subtract(l2,l1), 	
        //   math.subtract(
        //   	math.dotMultiply(l1,l3),
        //    math.dotMultiply(l2,l3))]];
        
        //abcd = math.multiply(
        //	math.inv(mat1),
        //  mat2);
         
         //var a 	= abcd[0][0];
         //var b 	= abcd[0][1];
         //var c 	= abcd[1][0];
         //var d 	= abcd[1][1];
        
         return [a,b,c,d]
      }
      
			// =============================================================================
      // FUNCTION CMULTIPLY
      // =============================================================================
      
      function cmultiply(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        
        var c1_len = c1.length;
        if (c1_len == undefined) {c1_len = 1;}
        
        var c2_len = c2.length;
        if (c2_len == undefined) {c2_len = 1;}
        
      	if (c1_len == 1 && c2_len == 1) { // Both variables are scalars
          var result 		= math.complex(
          	num1.re*num2.re-num1.im*num2.im, 
            num1.re*num2.im+num1.im*num2.re);
        } else if (c1_len > 1 && c2_len > 1) { // Both variables are vectors
        	var result 		= []
          for (i = 0; i < num1.length; i++) {
            result.push(math.complex(
              num1[i].re*num2[i].re-num1[i].im*num2[i].im, 
              num1[i].re*num2[i].im+num1[i].im*num2[i].re));
          }
        } else if (c1_len > 1 && c2_len == 1) { // The first variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num1.length; i++) {
            result.push(math.complex(
              num1[i].re*num2.re-num1[i].im*num2.im, 
              num1[i].re*num2.im+num1[i].im*num2.re));
          } 
        } else if (c1_len == 1 && c2_len > 1) { // The second variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num2.length; i++) {
            result.push(math.complex(
              num1.re*num2[i].re-num1.im*num2[i].im, 
              num1.re*num2[i].im+num1.im*num2[i].re));
          } 
        } else {
          var result 		= False;
          }

      return result;   
      }
      
			// =============================================================================
      // FUNCTION CDIVIDE
      // =============================================================================
      
      function cdivide(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        
        var c1_len = c1.length;
        if (c1_len == undefined) {c1_len = 1;}
        
        var c2_len = c2.length;
        if (c2_len == undefined) {c2_len = 1;}
        
      	if (c1_len == 1 && c2_len == 1) { // Both variables are scalars
          var denom 		= num2.im * num2.im + num2.re * num2.re;
          var real 			= (num1.re * num2.re + num1.im * num2.im) /denom;
          var imaginary = (num2.re * num1.im - num1.re * num2.im) /denom; 
          var result 		= math.complex(real, imaginary);
        } else if (c1_len > 1 && c2_len > 1) { // Both variables are vectors
        	var result 		= []
          for (i = 0; i < num1.length; i++) {
            var denom 		= num2[i].im * num2[i].im + num2[i].re * num2[i].re;
            var real 			= (num1[i].re * num2[i].re + num1[i].im * num2[i].im) /denom;
            var imaginary = (num2[i].re * num1[i].im - num1[i].re * num2[i].im) /denom; 
            result.push(math.complex(real, imaginary));
          }
        } else if (c1_len > 1 && c2_len == 1) { // The first variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num1.length; i++) {
            var denom 		= num2.im * num2.im + num2.re * num2.re;
            var real 			= (num1[i].re * num2.re + num1[i].im * num2.im) /denom;
            var imaginary = (num2.re * num1[i].im - num1[i].re * num2.im) /denom; 
            result.push(math.complex(real, imaginary));
          } 
        } else if (c1_len == 1 && c2_len > 1) { // The second variable is a vector
        	var result 		= [];
          var i 				= 0;
          for (i = 0; i < num2.length; i++) {
            var denom 		= num2[i].im * num2[i].im + num2[i].re * num2[i].re;
            var real 			= (num1.re * num2[i].re + num1.im * num2[i].im) /denom;
            var imaginary = (num2[i].re * num1.im - num1.re * num2[i].im) /denom; 
            result.push(math.complex(real, imaginary));
          } 
        } else {
          var result 		= False;
          }

      return result;   
      }
      
      function cadd(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        var real = (num1.re + num2.re);
        var imaginary = (num1.im + num2.im); 
      return math.complex(real, imaginary);   
      }
      
      function csub(c1, c2) {
        var num1, num2;
        num1 = c1;
        num2 = c2;
        var real = (num1.re - num2.re);
        var imaginary = (num1.im - num2.im); 
      return math.complex(real, imaginary);   
      }
      
      function cpow(c, exp) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push(math.pow(c[i],exp))
        }
      return res;   
      }
      
      function creal(c) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push(c[i].re)
        }
      return res;   
      }
      
      function cimag(c) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push(c[i].im)
        }
      return res;   
      }
      
      function rvec_to_cvec(c) {
        var res 	= [];
        for (i = 0; i < c.length; i++) {
          res.push([c[i]])
        }
      return res;   
      }
      
      
			// =============================================================================
      // UPDATE MODEL AND FIGURE
      // =============================================================================
      
      function update() {
        
      	p1 				= [
        	(d3.select("#point1").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#point1").attr("y")-height_offset)/height/0.9*2-1];
        p2 				= [
        	(d3.select("#point2").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#point2").attr("y")-height_offset)/height/0.9*2-1];
        w1 				= [
        	(d3.select("#well1").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#well1").attr("y")-height_offset)/height/0.9*2-1];
        w2 				= [
        	(d3.select("#well2").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#well2").attr("y")-height_offset)/height/0.9*2-1];
     
        z1 				= math.complex(
        	(d3.select("#zp1").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#zp1").attr("y")-height_offset)/height/0.9*2-1);
        z2 				= math.complex(
        	(d3.select("#zp2").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#zp2").attr("y")-height_offset)/height/0.9*2-1);
        z3 				= math.complex(
        	(d3.select("#zp3").attr("x")-width_offset)/width/0.9*2-1,
          (d3.select("#zp3").attr("y")-height_offset)/height/0.9*2-1);
        
        [a,b,c,d] = get_moebius_coefficients(z1,z2,z3,l1,l2,l3);
        
        p4 = moebius(z4,a,b,c,d);
        
       	d3.select("#zp4").attr("cx", (math.re(p4)+1)/2*width*0.9+width_offset);
        d3.select("#zp4").attr("cy", (math.im(p4)+1)/2*height*0.9+height_offset);
        
        vertices 	= linspace_complex(p1, p2, vres);
        zc 			 	= linspace_complex_segments(p1, p2, vres);
        nvec 			= normal_vector(p1,p2);
        A 				= noflowboundary_gradients(vertices,zc,nvec);
        solvec 		= solution_vector(zc,nvec);
        s 				= solve_linear_system(A,solvec);
        output 		= evaluate(grid, vertices, s);
        
        
        // Remove the previous elements
        d3.select("#line_element").remove();
        svg.selectAll("path").remove();


        // array of threshold values 
        
        //thresholds = d3.range(
        //d3.min(output) - (d3.max(output)- d3.min(output))/21,
        //d3.max(output) + (d3.max(output)- d3.min(output))/21*2,
        //(d3.max(output)- d3.min(output))/21);
        
        // Get min and max values inside the unit disk
        minval 	= +9999999
        maxval  = -9999999
        for (i = 0; i < output.length; i++) {
          if (math.abs(grid[i]) < 1) {
            if (output[i] < minval) {
              minval = output[i]
            }
            if (output[i] > maxval) {
              maxval = output[i]
            }
          }
        }
        
        // array of threshold values 
        thresholds = d3.range(
          minval*1.1,
          maxval,
          (maxval - minval)/21);

				color = d3.scaleLinear()
          .domain(d3.extent(thresholds))
          .interpolate(function() { return d3.interpolateRgbBasis(["#152d3b","#295c79","#4794c1","#c3e7f9"])});  

        // initialise contours
        contours = d3.contours()
            .size([y_res, x_res])
            .thresholds(thresholds)
            (output);

        // make and project the contours
        svg.selectAll("path")
            .data(contours)
            .enter().append("path")
                .attr("d", d3.geoPath(projection))
                .attr("fill", function(d) { return color(d.value); })
                .attr("clip-path", "url(#clipCircle)")

        svg
            .append("line")
            .attr("x1",(p1[0]+1)/2*width*0.9+width_offset)
            .attr("x2",(p2[0]+1)/2*width*0.9+width_offset)
            .attr("y1",(p1[1]+1)/2*height*0.9+height_offset)
            .attr("y2",(p2[1]+1)/2*height*0.9+height_offset)
            .attr("stroke", "#333333")
            .attr("stroke-width", 10)
            .attr("stroke-linecap","round")
            .attr("id","line_element");
        
        border.raise()
        point1.raise()
        point2.raise()
        well1.raise()
        well2.raise()
        
        zp1.raise()
        zp2.raise()
        zp3.raise()
        zp4.raise()
        
      return;   
      }
      

			// =============================================================================
      // PLOT THE CONTOURS
      // =============================================================================

      // set x and y scale to maintain 1:1 aspect ratio  
      // Extract the width and height that was computed by CSS.

      
      var scaling 	= math.min([width/x_res,height/y_res]);
      //print(scaling)
      
      var projection = d3.geoTransform({
          point: function(px, py) {
              this.stream.point(px*scaling*0.9+width_offset, py*scaling*0.9+height_offset);
          }
      });
      
			
        
      var svg = d3.select("#click") // This selects the div
          .attr("width", width) // This defines the canvas' width
          .attr("height", height) // This defines the canvas' height
      
      // Get min and max values inside the unit disk
      var minval 	= +9999999
      var maxval  = -9999999
      for (i = 0; i < output.length; i++) {
      	if (math.abs(grid[i]) < 1) {
        	if (output[i] < minval) {
          	minval = output[i]
          }
        	if (output[i] > maxval) {
          	maxval = output[i]
          }
        }
      }
      
      // array of threshold values 
      var thresholds = d3.range(
      	minval*1.1,
        maxval,
        (maxval - minval)/21);
      //var thresholds = d3.range(
      //	d3.min(output),
      //  d3.max(output),
      //  (d3.max(output)- d3.min(output))/21);

      // color scale  
      var color = d3.scaleLinear()
          .domain(d3.extent(thresholds))
          .interpolate(function() { return d3.interpolateRgbBasis(["#152d3b","#295c79","#4794c1","#c3e7f9"])});  
          //.interpolate(function() { return d3.interpolateRdBu; });  
          
      // initialise contours
      var contours = d3.contours()
          .size([y_res, x_res])
          .thresholds(thresholds)
          (output);


			//data = [];
      //for (i = 0; i < grid.length; i++) {
      //	data.push([math.re(grid[i])*40,math.im(grid[i])*40,output[i]])
      //}

			//var triangulation = d3.Delaunay.from(data /*, d => d.x, d => d.y */)

			//var contours = d3.tricontour()
      //    .thresholds(thresholds)
      //    .triangulate(() => triangulation)
      //    (data);


      // Resize the drag circles
      d3.select("#dragcircle").attr("r",20*height/450)
      
      // Get the clipping circle
      svg.append("clipPath")
        .attr("id", "clipCircle")
      .append("circle")
        .attr("r", width/2*0.9)
        .attr("cx", width/2*0.9+width_offset)
        .attr("cy", height/2*0.9+height_offset);
      
      // make and project the contours
      svg.selectAll("path")
          .data(contours)
          .enter().append("path")
              .attr("d", d3.geoPath(projection))
              .attr("fill", function(d) { return color(d.value); })
              .attr("clip-path", "url(#clipCircle)")
              
			svg
      		.append("line")
          .attr("x1",(p1[0]+1)/2*width*0.9+width_offset)
          .attr("x2",(p2[0]+1)/2*width*0.9+width_offset)
          .attr("y1",(p1[1]+1)/2*height*0.9+height_offset)
          .attr("y2",(p2[1]+1)/2*height*0.9+height_offset)
          .attr("stroke", "#333333")
          .attr("stroke-width", 10*height/450)
          .attr("stroke-linecap","round")
          .attr("id","line_element");
          
			var border = svg.append("circle")
        .attr("r", width/2*0.9)
        .attr("cx", width/2*0.9+width_offset)
        .attr("cy", height/2*0.9+height_offset)
        .attr("fill","None")
        .attr("stroke", "#666666")
        .attr("stroke-width", 3*height/450)
        .attr("id","border");;
          
			// =============================================================================
      // CREATE THE POINTER
      // =============================================================================
          
      var point1 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (p1[0]+1)/2*width*0.9+width_offset)
          .attr("y", (p1[1]+1)/2*height*0.9+height_offset)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","point1");
          
      var point2 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (p2[0]+1)/2*width*0.9+width_offset)
          .attr("y", (p2[1]+1)/2*height*0.9+height_offset)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","point2");
          
      var well1 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (w1[0]+1)/2*width*0.9+width_offset)
          .attr("y", (w1[1]+1)/2*height*0.9+height_offset)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","well1");
          
      var well2 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (w2[0]+1)/2*width*0.9+width_offset)
          .attr("y", (w2[1]+1)/2*height*0.9+height_offset)
          .attr("fill", "#c3e7f9")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","well2");
          
      var zp1 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (math.re(z1)+1)/2*width*0.9+width_offset)
          .attr("y", (math.im(z1)+1)/2*height*0.9+height_offset)
          .attr("lim_lo", 0)
          .attr("lim_hi", 1)
          .attr("flip", 0)
          .attr("fill", "#295c79")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","zp1");
          
      var zp2 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (math.re(z2)+1)/2*width*0.9+width_offset)
          .attr("y", (math.im(z2)+1)/2*height*0.9+height_offset)
          .attr("lim_lo", 0)
          .attr("lim_hi", 1)
          .attr("flip", 0)
          .attr("fill", "#295c79")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","zp2");
          
      var zp3 = svg.append("use")
          .attr("href", "#pointer")
          .attr("x", (math.re(z3)+1)/2*width*0.9+width_offset)
          .attr("y", (math.im(z3)+1)/2*height*0.9+height_offset)
          .attr("lim_lo", 0)
          .attr("lim_hi", 1)
          .attr("flip", 0)
          .attr("fill", "#295c79")
          .attr("stroke", "#333333")
          .attr("stroke-width", 5*height/450)
          .attr("id","zp3");
          
			var zp4 = svg.append("circle")
        .attr("r", 10*height/450)
        .attr("cx", (math.re(p4)+1)/2*width*0.9+width_offset)
        .attr("cy", (math.im(p4)+1)/2*height*0.9+height_offset)
        .attr("fill", "#295c79")
        .attr("stroke", "#666666")
        .attr("opacity", 0.5)
        .attr("stroke-width", 2.5*height/450)
        .attr("id","zp4");
        
          //.attr("transform","scale("+string(0.5*height/450)+")");

      var deltaX, deltaY;

      var dragHandler = d3.drag()
          .on("start", function () {
          		dragging = true;
              var current = d3.select(this);
              deltaX = current.attr("x") - d3.event.x;
              deltaY = current.attr("y") - d3.event.y;
              
              // If this element is an edge point, get the difference to its neighbours
              if (current.attr("id") == "zp1" ||
              		current.attr("id") == "zp2" ||
                  current.attr("id") == "zp3") {
              
              	// Find the neighbours
                if (current.attr("id") == "zp1") {
                  neighbour1 	= d3.select("#zp2")
                  neighbour2 	= d3.select("#zp3")
                  }
                if (current.attr("id") == "zp2") {
                  neighbour1 	= d3.select("#zp1")
                  neighbour2 	= d3.select("#zp3")
                  }
                if (current.attr("id") == "zp3") {
                  neighbour1 	= d3.select("#zp1")
                  neighbour2 	= d3.select("#zp2")
                  }
                  
                  // Get the angles
                  ac 					= math.atan2(current.attr("y")-height/2,current.attr("x")-width/2)
                  anb1 				= math.atan2(neighbour1.attr("y")-height/2,neighbour1.attr("x")-width/2)
                  anb2 				= math.atan2(neighbour2.attr("y")-height/2,neighbour2.attr("x")-width/2)
                  
                  // Find out if the the borders need flipping or not
                  if (((ac > anb1) && (ac > anb2)) || ((ac < anb1) && (ac < anb2))) {
                  	
                    // The current point is either larger or lower than both neighbours
                    flip 			= 1
                    
                    // Find the upper and lower limit
                    lim_lo 		= math.max([anb1,anb2]) + 0.5
                    lim_hi 		= math.min([anb1,anb2]) - 0.5
                    
                  } else {
                  
                  	// the current point is between the neighbours
                    flip 			= 0
                    
                    // Find the upper and lower limit
                    lim_lo 		= math.min([anb1,anb2]) + 0.5
                    lim_hi 		= math.max([anb1,anb2]) - 0.5
                    
                  }
                  
                  // Wrap around if outside -pi to +pi
                  if (lim_lo < -math.PI) {
                  	lim_lo 	+= math.PI
                  }
                  if (lim_hi > math.PI) {
                  	lim_hi	-= math.PI
                  }
               
               		// Save the results to the variable
                  d3.select(this)["_groups"][0][0].flip 	= flip
                  d3.select(this)["_groups"][0][0].lim_lo	= lim_lo
                  d3.select(this)["_groups"][0][0].lim_hi	= lim_hi
                  
               		//current.attr("flip") 		= flip
                  //current.attr("lim_lo") 	= lim_lo
                  //current.attr("lim_hi") 	= lim_hi
   
              }
              
              
          })
          .on("drag", function () {
          		var movex = d3.event.x + deltaX;
              var movey = d3.event.y + deltaY;
              var current = d3.select(this);
  
              
              // Get the nopixel coordinates
              movex_np 	= (movex-width_offset)/width/0.9*2-1
              movey_np 	= (movey-height_offset)/height/0.9*2-1
              
							// Get the norm of the move
              norm 	= math.sqrt(
               	math.pow(movex_np,2) + math.pow(movey_np,2))
              
              if (current.attr("id") == "zp1" ||
              		current.attr("id") == "zp2" ||
                  current.attr("id") == "zp3") {
                  
                // Convert it to an angle
                alpha = math.atan2(movey_np,movex_np)
                
                // Extract the variables from earlier
                flip 		= d3.select(this)["_groups"][0][0].flip
                lim_lo 	= d3.select(this)["_groups"][0][0].lim_lo
                lim_hi 	= d3.select(this)["_groups"][0][0].lim_hi
                
                // Check for boundaries
                if (flip == 0) {
                
                  if ((alpha < lim_lo) || (alpha > lim_hi)) {

                    // Correct the angle
                    if (math.sqrt(
                          math.add(
                            math.pow(
                              math.cos(alpha)-math.cos(lim_lo),
                              2),
                            math.pow(
                              math.sin(alpha)-math.sin(lim_lo),
                              2))) < math.sqrt(
                          math.add(
                            math.pow(
                              math.cos(alpha)-math.cos(lim_hi),
                              2),
                            math.pow(
                              math.sin(alpha)-math.sin(lim_hi),
                              2)))) {

                     // Point is closer to lim_lo than lim_hi
                     alpha 		= lim_lo

                     } else {

                     // Point is closer to lim_hi than lim_lo
                     alpha 		= lim_hi

                     }

                  }
                  
                } else {
                
                  if ((alpha < lim_lo) && (alpha > lim_hi)) {

                    // Correct the angle
                    if (math.sqrt(
                          math.add(
                            math.pow(
                              math.cos(alpha)-math.cos(lim_lo),
                              2),
                            math.pow(
                              math.sin(alpha)-math.sin(lim_lo),
                              2))) < math.sqrt(
                          math.add(
                            math.pow(
                              math.cos(alpha)-math.cos(lim_hi),
                              2),
                            math.pow(
                              math.sin(alpha)-math.sin(lim_hi),
                              2)))) {

                     // Point is closer to lim_lo than lim_hi
                     alpha 		= lim_lo

                     } else {

                     // Point is closer to lim_hi than lim_lo
                     alpha 		= lim_hi

                     }

                  }
                
                }
                
                
                movex_np = math.cos(alpha)
                movey_np = math.sin(alpha)

                
                
                } else {
                
                if (norm > 1) {
                  movex_np 	= movex_np/norm
                  movey_np 	= movey_np/norm
                }
              
              }
              
 
              
              // Convert this back to pixels
              movex 	= (movex_np+1)/2*width*0.9+width_offset
							movey 	= (movey_np+1)/2*height*0.9+height_offset
              
              current
                  .attr("x", movex)
                  .attr("y", movey);
                  
              update();
                
                  
          })
          .on("end", function () {
              update();
              dragging = false;
          });


      dragHandler(svg.selectAll("use"));




			
      // ===========================================================================
      // This function shifts the contours
      // ===========================================================================
      
      var increment = 0;
      
      function shift_contours() {


				threshdif = thresholds[1]-thresholds[0];
				increment += threshdif/10
        
        
        if (increment >= threshdif){
        	increment = 0;
        }

        // initialise contours
        contours = d3.contours()
            .size([y_res, x_res])
            .thresholds(math.subtract(thresholds,increment))
            (output);

  			//contours = d3.tricontour()
        //		.triangulate(() => triangulation)
        //    .thresholds(math.subtract(thresholds,increment))
        //    (data)
        
				svg.selectAll("path").remove();
        
        // make and project the contours
        svg.selectAll("path")
            .data(contours)
            .enter().append("path")
                .attr("d", d3.geoPath(projection))
                .attr("fill", function(d) { return color(d.value); })
                .attr("clip-path", "url(#clipCircle)")

				d3.select("#line_element").raise()
        border.raise()
        point1.raise()
        point2.raise()
				well1.raise()
        well2.raise()
				zp1.raise()
        zp2.raise()
        zp3.raise()
        zp4.raise()
				return }


  		// Make the contour flow
      var repater = setInterval(function() {
      	if (dragging == false) {shift_contours();}
      }, 25);
      
  
  
  
    </script>
  </body>

</html>









## Elements

<img align="left" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/01_uniform.png" width="15%">

### Uniform base flow
To provide a background potential or unidirectional regional flow, the simplest option is to use uniform flow, specified by the AEM toolbox' object `ElementUniformBase`. This element requires the specification of a direction in radians, a minimum and maximum hydraulic head, and a background hydraulic conductivity.
<br /><br />

<img align="right" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/02_moebius.png" width="15%">

### Möbius base flow
Möbius base flow provides a way to implement more intricate regional flow, allowing for curvature, divergence, and convergence. This type of flow is specified by the AEM toolbox' object `ElementMoebiusBase`. This element requires the specification of three control points' direction in radians, a minimum and maximum hydraulic head, and a background hydraulic conductivity.
<br /><br />

<img align="left" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/03_extraction_well.png" width="15%">

### Extraction or injection well
Injection or extraction wells - or any other type of pointwise flow - can be implemented using the `ElementWell` object. This element requires the specification of a position, a positive or negative strength value, and a well radius. Alternatively to a strength value, this element can also adjust its strength to induce a desired drawdown on the flow field.
<br /><br />

<img align="right" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/04_inhomogeneity.png" width="15%">

### Polygonal inhomogeneity
Zonal inhomogeneities in the aquifer's hydraulic conductivity can be represented using the `ElementInhomogeneity` object. This element requires the specification of a hydraulic conductivity value, as well as a closed or open polygon defining its extent.
<br /><br />

<img align="left" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/05_fixed_head_boundary.png" width="15%">

### Prescribed head boundary / River
Prescribed head boundary conditions or rivers can be created using the `ElementHeadBoundary` object. This line element enforces a specified hydraulic head along its path. It requires the specification of its vertices and corresponding head values. It also allows for the implementation of a uniform or spatially varying connectivity value which can limit its influence on the water table.
<br /><br />

<img align="right" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/06_no_flow_boundary.png" width="15%">

### No-flow boundary
No-flow boundaries from sheet pile walls or impermeable formations can be created using the `ElementNoFlowBoundary` object. This line element requires only the specification of the vertices along its path and can be either closed or open.
<br /><br />

<img align="left" src="https://raw.githubusercontent.com/MaxRamgraber/Simple-AEM-Toolbox/main/images/07_area_sink.png" width="15%">

### Area sources or sinks
Areal recharge or water extraction can be represented using the `ElementAreaSink` object. This element adds or removes water according to a specified range inside its polygon. It requires the specification of its polygons and a positive or negative strength value. Not that the stream function component of the complex potential is not valid inside an area source or sink.
<br /><br />

## Troubleshooting

**Q: The model seems to create singularities, predicting very high or low water tables at certain isolated locations. What did I do wrong?**

A: This usually happens if the model attempts to evaluate the complex potential Ω directly on an element. This can happen because because two elements share a line segment or because one of the evaluation points lies on an a line segment. Make sure that the elements do not share direct borders, for example by offsetting them by a minuscule amount (e.g., 1E-10). I have implemented protections against this for some but not all elements: inhomogeneity elements, for example, are automatically shrunk by a negligible amount. Also, you should make sure that no inhomogeneities or no-flow boundaries intersect.

**Q: There are still strange artefacts along my no-flow boundaries or inhomogeneities. What happened?**

A: If you have tried the solutions in the answer above and the issue persists, try increasing the resolution of the element by increasing the element's `segments`. Most of the constant-strength line elements I used here require sufficient resolution to induce the desired effect. It is difficult to predict how large this resolution should be in advance.

**Q: Why are there discontinuities in the stream function away from my wells, head boundaries, or area sinks?**

A: These are so-called branch cuts, discontinuities in the stream function (the imaginary part of the complex potential) which arise whenever an element adds or removes water from the system. An example can be seen in the last image created in the [basic tutorial](https://github.com/MaxRamgraber/Simple-AEM-Toolbox/tree/main/Tutorials/Tutorial%2001%20Basic%20AEM). Unfortunately, there is no way to avoid them. For particle tracking purposes (when using the stream function, not the hydraulic potential), I advise capitalizing on the fact that these branch cuts occur in predictable patterns and explicitly account for them in your tracking routine. For plotting purposes, I recommend not using ready-made contouring functions, but instead plotting the streamlines with with a particle tracking routine.

**Q: The stream function appears discontinuous between the inside and outside of my area source or sink. Why?**

A: The stream function is not valid inside of area sinks or sources. If you want to plot the stream function anyways, I recommend masking the stream function inside the area sink's or source's polygon, similar to what I did in the example area sink under the 'Elements' heading above.
