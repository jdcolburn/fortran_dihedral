module data_types

	! declaration of geometric objects
	
	type :: point
		real, 			dimension(3) :: coords			 !cartesian coordinates
	end type point                                       
	                                                     
	type :: vector                                       
		character, 		dimension(2) :: name			 !names of two component points
		type( point  ), dimension(2) :: component_point	 !two component points
		real						 :: length			 !scalar length
	end type vector                                      
	                                                     
	type :: plane                                        
		character, 		dimension(3) :: name			 !names of three component points
		type( point  ), dimension(3) :: component_point	 !three component points
		type( vector ), dimension(2) :: component_vector !two component vectors
		type( vector ) 				 :: normal			 !normal vector
	end type plane                   
	                                 
	! declaration of chemical objects  
	
	type :: atom                     
		character ( len=1 ) 		 :: label
		integer 					 :: id
		type( point )				 :: position
	end type atom
	
	type :: torsion
		character, 		dimension(4) :: name
		type( plane ), 	dimension(2) :: component_plane
		real						 :: angle
	end type torsion	
	
end module data_types	

module functions

	use data_types
	implicit none
	
	contains
	
	function pythagoras( a, b )
	! returns hypotenuse - root square of component coords 
	
		TYPE( point )	:: a, b			! input 
		real			:: pythagoras 	! output
		
			pythagoras = sqrt ( &
			(( a%coords(1) - b%coords(1) )**2) + &
			(( a%coords(2) - b%coords(2) )**2) + &
			(( a%coords(3) - b%coords(3) )**2) &
			)
		
	end function pythagoras
	
	function cross_product( u, v )
	! returns cross product of two vectors (also coefficients for a vector perpendicular to a plane)
	
		type( vector )	:: u, v
		type( point  )  :: cross_product
		
			cross_product%coords(1) = ( u%component_point(2)%coords(2) * v%component_point(2)%coords(3) - &
										u%component_point(2)%coords(3) * v%component_point(2)%coords(2) )			
			cross_product%coords(2) = ( u%component_point(2)%coords(3) * v%component_point(2)%coords(1) - &
										u%component_point(2)%coords(1) * v%component_point(2)%coords(3) )
			cross_product%coords(3) = ( u%component_point(2)%coords(1) * v%component_point(2)%coords(2) - &
										u%component_point(2)%coords(2) * v%component_point(2)%coords(1) )

	end function cross_product
	
	function define_vector( i, j )
	! creates and defines all attributes for a vector datatype
		
		type( atom )	:: i, j			 ! input 
		type( vector )	:: define_vector ! output
		
			!assign name
			define_vector%name(1) = ( i%label )
			define_vector%name(2) = ( j%label )
			
			!assign component points centered at origin
			define_vector%component_point(1)%coords(1) = 0
			define_vector%component_point(1)%coords(2) = 0
			define_vector%component_point(1)%coords(3) = 0		
			
			define_vector%component_point(2)%coords(1) = ( i%position%coords(1) - j%position%coords(1) )
			define_vector%component_point(2)%coords(2) = ( i%position%coords(2) - j%position%coords(2) )
			define_vector%component_point(2)%coords(3) = ( i%position%coords(3) - j%position%coords(3) )
		
			!assign length
			define_vector%length = pythagoras( define_vector%component_point(1), define_vector%component_point(2) )
	
	end function define_vector
	
	function define_plane( i, j, k )
	! creates and defines all attributes for a plane datatype
		
		type( atom )  :: i, j, k       ! input 
		type( plane ) :: define_plane  ! output
		
			!assign name
			define_plane%name(1) = ( i%label )
			define_plane%name(2) = ( j%label )
			define_plane%name(3) = ( k%label )
			
			!assign three component points
			define_plane%component_point(1) = ( i%position )
			define_plane%component_point(2) = ( j%position )
			define_plane%component_point(3) = ( k%position )
			
			!assign two component vectors
			define_plane%component_vector(1) = define_vector( i, j )
			define_plane%component_vector(2) = define_vector( i, k )
			
			!assign normal vector
			define_plane%normal%component_point(1)%coords(1) = 0
			define_plane%normal%component_point(1)%coords(2) = 0
			define_plane%normal%component_point(1)%coords(3) = 0		
			define_plane%normal%component_point(2) = cross_product( define_plane%component_vector(1), define_plane%component_vector(2) )
			
	end function define_plane	

	function get_dihedral( i, j, k, l )
	! calculates the dihedral angle between four atoms
	
		type( atom ) 	 :: i, j, k, l    ! input 
		type( torsion )  :: get_dihedral  ! output
	
			!assign name
			get_dihedral%name(1) = ( i%label )
			get_dihedral%name(2) = ( j%label )
			get_dihedral%name(3) = ( k%label )
			get_dihedral%name(4) = ( l%label )
			
			!assign two component planes
			get_dihedral%component_plane(1) = define_plane( i, j, k )
			get_dihedral%component_plane(2) = define_plane( j, k, l )
			
			!calculate angle
			get_dihedral%angle = acos ( ( &
			(get_dihedral%component_plane(1)%normal%component_point(2)%coords(1) * &
			 get_dihedral%component_plane(2)%normal%component_point(2)%coords(1)) + &								 
			(get_dihedral%component_plane(1)%normal%component_point(2)%coords(2) * &
			 get_dihedral%component_plane(2)%normal%component_point(2)%coords(2)) + &			
			(get_dihedral%component_plane(1)%normal%component_point(2)%coords(3) * &
			 get_dihedral%component_plane(2)%normal%component_point(2)%coords(3)) ) &			
			/ ( &		
			pythagoras( get_dihedral%component_plane(1)%normal%component_point(1), &
			get_dihedral%component_plane(1)%normal%component_point(2)) * &
			pythagoras( get_dihedral%component_plane(2)%normal%component_point(1), &
			get_dihedral%component_plane(2)%normal%component_point(2)) &
			))
								 
	end function get_dihedral
	
end module functions

program dihedral
	
	use data_types
	use functions
	implicit none
 
 !variabe declaration
 integer :: i, j, k, l
 type( vector ) :: define_vector
 type( plane )  :: define_plane
 type( torsion ) :: torsion_ijkl
 type( atom ), dimension(8) :: molecule 
 
	molecule(1)%id=1
	molecule(1)%label="C"
	molecule(1)%position%coords(1)= 0.762209
	molecule(1)%position%coords(2)= 0.000000
	molecule(1)%position%coords(3)= 0.000000
	
	molecule(2)%id=2
	molecule(2)%label="C"
	molecule(2)%position%coords(1)= 0.000000
	molecule(2)%position%coords(2)= 0.000000
	molecule(2)%position%coords(3)=-0.762209
	
	molecule(3)%id=3
	molecule(3)%label="H"
	molecule(3)%position%coords(1)= 0.000000
	molecule(3)%position%coords(2)= 1.018957 
	molecule(3)%position%coords(3)= 1.157229
	
	molecule(4)%id=4
	molecule(4)%label="H"
	molecule(4)%position%coords(1)=-0.882443
	molecule(4)%position%coords(2)=-0.509479
	molecule(4)%position%coords(3)= 1.157229
	
	molecule(5)%id=5
	molecule(5)%label="H"
	molecule(5)%position%coords(1)= 0.882443
	molecule(5)%position%coords(2)=-0.509479
	molecule(5)%position%coords(3)= 1.157229
	
	molecule(6)%id=6
	molecule(6)%label="H"
	molecule(6)%position%coords(1)= 0.000000
	molecule(6)%position%coords(2)=-1.018957 
	molecule(6)%position%coords(3)=-1.157229
	
	molecule(7)%id=7
	molecule(7)%label="H"
	molecule(7)%position%coords(1)=-0.882443
	molecule(7)%position%coords(2)= 0.509479
	molecule(7)%position%coords(3)=-1.157229
	
	molecule(8)%id=8
	molecule(8)%label="H"
	molecule(8)%position%coords(1)= 0.882443
	molecule(8)%position%coords(2)= 0.509479
	molecule(8)%position%coords(3)=-1.157229

	
print*, molecule(1)%id, molecule(1)%label
print*, molecule(2)%id, molecule(2)%label
print*, molecule(3)%id, molecule(3)%label
print*, molecule(4)%id, molecule(4)%label
print*, molecule(5)%id, molecule(5)%label
print*, molecule(6)%id, molecule(6)%label
print*, molecule(7)%id, molecule(7)%label
print*, molecule(8)%id, molecule(8)%label

print*, "Select four atom IDs: "	
read*, i
read*, j
read*, k
read*, l

torsion_ijkl = get_dihedral( molecule(i), molecule(j), molecule(k), molecule(l) )

print*, "Torsional angle for ", torsion_ijkl%name, " is: ", torsion_ijkl%angle, "degrees."

end program dihedral
