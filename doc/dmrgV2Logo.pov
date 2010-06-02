#include "colors.inc"
#include "finish.inc"
#include "textures.inc"
background { color White }

camera {
	location <-2,-15,-1>
	look_at <-2,0,0>
	rotate 10*x
	rotate <0,0,clock*360>
} 

light_source {  <2,8,-4> spotlight point_at<2,0,0> }
light_source { <-2,8,-4> White }
light_source {  <4,0,-2>  }
light_source { <-4,0,-2> White }
plane { z, 0
	//rotate -15*x
	translate 1.5*z
	translate 30*y
	texture {
		Water
	}
}

plane { z, 0
	//rotate 15*x
	translate 30*y
	translate -5*z
	texture {Bright_Blue_Sky
		//pigment {color SkyBlue}
		finish{Luminous}
	}
}

box {
	<0, 0, 0>
	<0,200,50>
	translate -100*y
	translate -25*z

	texture { 
		pigment {Gray}
		finish {Mirror} 
	}

}
#declare MyDifference = 
difference {
	box {
		<0, 0, 0>
		<0,200,48>
     		texture { DMFWood6 }
	}
	box {
		<0, 0, 0>
		<0,200,48>
      		texture { DMFWood6 }
  	}
}

cylinder
{
	<0,0,0>, <12,0,0>, 0.5
	translate 11*y
	translate -15*x
	texture { 
   		pigment {color Red} 
   		finish {Dull} 
	}
}

sphere {
	<0, 0, 0>, 2
	translate 10*y
    
	translate -3*x
	texture { 
		pigment {color Red} 
		finish {Metal} 
	}

}

box {
	<0,0,0>
	<13,4,4>
	
	translate -20*x
	translate 10*y
	translate -2*z
	texture { 
		pigment {color Red} 
		finish {Dull} 
	}
} 
  
text {
	ttf "timrom.ttf" "DMRG++" 0.1, 0
	texture { 
		pigment {color Blue} 
		finish {Shiny }
	}

	rotate -x*80
	scale 4
	translate -14.7*x
	translate +10*y
	translate -2*z
}
   
