!!!=========================************PERSONAL PROGRAMMING PROJECT****************===========================================!!!
!!!=========CODE AUTHOR: PADMANABHA PAVAN CHANDRA VUNDURTHY(Mat.Nr.62750)(MSc. COMPUTATIONAL MATERIALS SCIENCE)========================!!!
!!!============================UNDER THE SUPERVISION OF DR. ARUN PRAKASH========================!!!
!!!============TITLE: BIG DATA VISUALIZATION OF LARGE AND ULTRA LARGE SCALE ATOMISTIC SIMULATION (UPTO 20 MILLION ATOMS) =================!!!
!!!=====FILES NECESSARY AS INPUT: PARAMETER FILE, SORTED DATA FILES CONTAINING THE PROPERTIES NECESSARY FOR ANALYSIS THAT ARE MENTIONED IN PARAMETER FILE====!!!  
!!!======================OUTPUT FILES: CONTAINS THE AVERAGE PROPERTIES OF THE GRAINS WITH DIFFERENT STRUCTURE TYPES=======================!!!
!!!===================THE FOLLOWING PROGRAM WORKS FOR 'N' NUMBER OF FILES AND ANY AMOUNT OF DATA PROVIDED FOR ANALYSIS==============!!!
PROGRAM Datai_avg
IMPLICIT NONE
	integer::error, num_files,a, num_props, b, i,j, errro,k,l,m
	character(120), dimension(:), allocatable:: filename_arr, prop_name_arr
	integer,dimension(:,:), allocatable::file_attributes_arr, prop_attributes_arr
	integer::c, atomid_filenum, d=1, e, f, g, h, span, counter, grnmax, structmax
	integer::n,o, p,q,r,s,t,u
	real, dimension(:,:), allocatable::iprop_arr, temp_arr
	real,dimension(:), allocatable::atomid_arr, vm_s
	real, dimension(:,:,:), allocatable::avg_props
	real::vm_v,x,y,z,xy,yz,zx,Macro_VM,summer=0.0
	!!================================!!
	!!*ABOVE DECLARED VARIABLES*!!
	!*Data from parameter File:*!
		!filename_arr->array contains the filenames, prop_name_arr->array contains property names
		!file_attributes_arr->array contains file numbers w.r.t property names
		!prop_attributes_arr->array contains the location of a property  
	!temp_arr-> temporary array of a predefined size that is determined from the parameter file 
	!iprop_arr-> the final array which collects the required data from all the files based on user's interests
	!avg_props->a 3D array containing the average properties for all the grains w.r.t structure types represented as a block
	!mACRO_VM -> Von Mises Stress
	!!================================!!	
	!***OPEN AND READ THE PARAMETER FILE FOR THE USER INPUT***!
	open (unit=201, file='Parameter_File_BigData.txt', iostat = error)
	if (error .ne. 0) then
            write(*,*) "file cannot be opened"
    end if
    read(201, *)
    read(201,*)num_files
	print*,"Number of Files==",num_files
    read(201, *)
    allocate(filename_arr(num_files))
    allocate(file_attributes_arr(3,num_files))
    !READ FOR STORING THE FILE NAMES AND FILE ATTRIBUTES
    a = 0
    do while(a<num_files)
		read(201,*)filename_arr(a+1), file_attributes_arr(:,a+1) 
		a = a+1
	end do
	a = 1
	do a = 1,num_files
		print*,"filename::==",filename_arr(a)
		print*, "file attributes::==",file_attributes_arr(:,a)
	end do
	read(201,*)
	read(201,*)
	read(201,*)num_props
	!!ALLOCATE THE ALLOCATABLE ARRAY DIMENSIONS FROM THE INPUT READ FROM PARAMETER FILE 
	allocate(prop_name_arr(num_props))	
	allocate(prop_attributes_arr(3,num_props))	
	b = 0
	do while(b<num_props)
		read(201,*)prop_name_arr(b+1), prop_attributes_arr(:,b+1)
		b = b+1
	end do
	close(201)
	!CLOSE THE PARAMETER FILE
	!!=========================================!!
	!!=========================================!!
	! ASSIGN ATOM ID AND CREATE AN ARRAY OF INTERESTED QUANTITIES
	!useful arrays: filename_arr, file_attribiutes_arr,prop_name_arr,prop_attribites_arr, num_props, num_files
	c = 1
	do while (c<=num_props)
		if (trim(prop_name_arr(c))=='atomID') then
			atomid_filenum = prop_attributes_arr(2,c)
			!IDENTIFYING THE MAXIMUM VALUE OF THE NUMBER OF COLUMNS & NUMBER OF ROWS IN THE FILES
			i = maxval(file_attributes_arr(2,:))		
			j = maxval(file_attributes_arr(1,:))-minval(file_attributes_arr(3,:))	
			allocate(temp_arr(i,j))		!!ALLOCATE A SIZE FOR TEMPORARY ARRAY W.R.T MAX. DIMENSIONS POSSIBLE
			!!INITIALISING THE ARRAY TO '0'!!
			l=1
			do while(l<=j)
				temp_arr(:,l)=0
				l=l+1
			end do
			!ALLOCATE ATOM_ID ARRAY DIMENSIONS
			allocate(atomid_arr(file_attributes_arr(1,atomid_filenum)-file_attributes_arr(3,atomid_filenum)))
			!ALLOCATE THE DIMENSIONS FOR IPROP ARRAY BASED ON THE NUMBER OF INTERESTED QUANTITIES FROM THE USER
			allocate(iprop_arr(13,size(atomid_arr)))
			!*OPEN AND READ THE FILE CONTAINING ATOM_ID *!
			open(unit=202,file = filename_arr(atomid_filenum), iostat = errro)	
			if (errro .ne. 0) then
					write(*,*) "file cannot be opened"
			end if	
			!!**READING THE FILE INTO A TEMPORARY ARRAY**!!
			k = 1
			do
				if (k<=file_attributes_arr(3,atomid_filenum))then
					read (202,*,end=10)
				else if (k<=file_attributes_arr(1,atomid_filenum)) then
					read (202,*,end=10) temp_arr(1:file_attributes_arr(2,atomid_filenum),k-file_attributes_arr(3,atomid_filenum))					
				end if
				k=k+1
			end do
10 continue
			close(202)
			!*CLOSED THE FILE*!
			!print*,temp_arr(:,3)
			!EXTRACTING THE ATOM ID ARRAY FROM TEMPORARY ARRAY
			atomid_arr = temp_arr(1,1:file_attributes_arr(1,atomid_filenum)-file_attributes_arr(3,atomid_filenum))		
			iprop_arr(1,:)=atomid_arr		!ASSIGNING ATOM_ID ARRAY TO THE FIRST COLUMN OF IPROP ARRAY
			!!!***START CHECKING FOR PARAMETER TAGS***!!!
			do while (d<=num_props)
				if (prop_attributes_arr(2,d)==atomid_filenum)then					
				!**-------------------------------------------**!
					if (trim(prop_name_arr(d))=='Coord') then
						iprop_arr(2:4,:)=temp_arr(prop_attributes_arr(3,d):prop_attributes_arr(3,d)+2,:)
						!print*,"=========", iprop_arr(1:4,:)
					else if (trim(prop_name_arr(d))=='struct_type') then
						iprop_arr(5,:)=temp_arr(prop_attributes_arr(3,d),:)
						!print*,"=========", iprop_arr(1:5,:)
					else if (trim(prop_name_arr(d))=='GrNum') then
						iprop_arr(6,:)=temp_arr(prop_attributes_arr(3,d),:)
						!print*,"=========", iprop_arr(1:6,:)
					else if (trim(prop_name_arr(d))=='Stress') then
						iprop_arr(7:12,:)=temp_arr(prop_attributes_arr(3,d):prop_attributes_arr(3,d)+5,:)
						!print*,"=========", iprop_arr(1:12,:)					
					else if (trim(prop_name_arr(d))=='Volume') then
						iprop_arr(13,:)=temp_arr(prop_attributes_arr(3,d),:)						
						!print*,"=========", iprop_arr(1:13,:)
					end if
				!**-------------------------------------------**!
				end if
				d = d+1
			end do
			!!!***END CHECKING FOR PARAMETER TAGS***!!!
		end if
		c = c+1
	end do
	!!===============================================================!!
	!!===============================================================!!
	!CHECK FOR ATOM ID FILE, AS THE TASKS TO BE PERFORMED IN THAT FILE ARE...
	!...ALREADY ACCOMPLISHED BY EXTRACTING OTHER PROPERTIES IN THAT FILE AND STORING IN IPROP ARRAY.
	!NOW, EXTRACT NECESSARY DATA FROM OTHER FILES IGNORING THE ATOM ID FILE
	e = 1
	do while (e<=num_files)		!!LOOP OVER THE NUMBER OF FILES 
		! IGNORE ATOM ID FILE
		if (e==atomid_filenum)then
			goto 25
		end if
		!!INITIALISE THE TEMP ARRAY TO '0' FOR AVOIDING ANY GARBAGE VALUES
		l=1
		do while(l<=j)
			temp_arr(:,l)=0
			l=l+1
		end do
		span = file_attributes_arr(1,e)-file_attributes_arr(3,e)	!EXTRACTS NUMBER OF ROWS
		!*OPEN THE FILE*!
		open(unit=203,file = filename_arr(e), iostat = errro)	
		if (errro .ne. 0) then
				write(*,*) "file cannot be opened"
		end if	
		!RE WRITE THE TEMP ARRAY WITH DATA FROM NEW FILE
		k = 0
		do
			k=k+1
			if (k<=file_attributes_arr(3,e))then
				read (203,*,end=20)
				!print*,n
			else if (k<=file_attributes_arr(1,e)) then
				read (203,*,end=20) temp_arr(1:file_attributes_arr(2,e),k-file_attributes_arr(3,e))
				
			end if
		end do
20 continue		
		close(203)
		!*CLOSED THE FILE*!
		!!=========================================================!!
		f = 1
		!LOOP OVER NUMBER OF PROPERTIES TO BE EXTRACTED
		do while(f<=num_props)			
			if (prop_attributes_arr(2,f)==e)then
				!!**CHECK FOR PARAMETER TAGS**!!
				!*CHECK FOR COORDINATES OF THE ATOM*!
				if (trim(prop_name_arr(f))=='Coord') then
					g = 1
					counter = 1
					do while(g<=size(atomid_arr))	
						h = counter
						do while(h<=span)							
							if (atomid_arr(g)==temp_arr(1,h))then
								iprop_arr(2:4,g)=temp_arr(prop_attributes_arr(3,f):prop_attributes_arr(3,f)+2,h)
								counter = h !COUNTS THE NO. OF LINES READ AND STARTS FROM THIS LINE IN THE NEXT ITERATION, SAVING TIME
								goto 35
							end if		
							h = h+1					
						end do
35 continue
						g = g+1
					end do
				!*CHECK FOR STRUCTURE TYPE*!
				else if (trim(prop_name_arr(f))=='struct_type') then
					g = 1
					counter = 1
					do while(g<=size(atomid_arr))
						h = counter
						do while(h<=span)
							if (atomid_arr(g)==temp_arr(1,h))then
								iprop_arr(5,g)=temp_arr(prop_attributes_arr(3,f),h)
								counter = h !COUNTS THE NO. OF LINES READ AND STARTS FROM THIS LINE IN THE NEXT ITERATION, SAVING TIME
								goto 45
							end if	
							h = h+1						
						end do
45 continue
						g = g+1
					end do
				!*CHECK FOR GRAIN NUMBER*!
				else if (trim(prop_name_arr(f))=='GrNum') then
					g = 1
					counter = 1
					do while(g<=size(atomid_arr))
						h = counter
						do while(h<=span)							
							if (atomid_arr(g)==temp_arr(1,h))then
								iprop_arr(6,g)=temp_arr(prop_attributes_arr(3,f),h)
								counter = h !COUNTS THE NO. OF LINES READ AND STARTS FROM THIS LINE IN THE NEXT ITERATION, SAVING TIME
								goto 55
							end if
							h = h + 1							
						end do
55 continue
						g = g+1
					end do
				!*CHECK FOR STRESS*!
				else if (trim(prop_name_arr(f))=='Stress') then
					g = 1
					counter = 1
					do while(g<=size(atomid_arr))	
						h = counter
						do while(h<=span)
							if (atomid_arr(g)==temp_arr(1,h))then
								iprop_arr(7:12,g)=temp_arr(prop_attributes_arr(3,f):prop_attributes_arr(3,f)+5,h)
								counter = h !COUNTS THE NO. OF LINES READ AND STARTS FROM THIS LINE IN THE NEXT ITERATION, SAVING TIME
								goto 65
							end if		
							h = h+1					
						end do
65 continue
						g = g+1
					end do
				!*CHECK FOR VOLUME*!				
				else if (trim(prop_name_arr(f))=='Volume') then
					g = 1
					counter = 1
					do while(g<=size(atomid_arr))
						h = counter
						do while(h<=span)
							if (atomid_arr(g)==temp_arr(1,h))then
								iprop_arr(13,g)=temp_arr(prop_attributes_arr(3,f),h)
								counter = h !COUNTS THE NO. OF LINES READ AND STARTS FROM THIS LINE IN THE NEXT ITERATION, SAVING TIME
								goto 75
							end if	
							h = h + 1						
						end do
75 continue
						g = g+1
					end do				
				end if
			end if				
			f = f+1
		end do
		!!============================================================!!		
25 continue
		e = e+1	
	end do
	!LOOPING OVER ALL THE FILES IS ACCOMPLISHED AND THE DATA NECESSARY FOR ANALYSIS IS EXTRACTED FROM FILES BASED ON USER'S INTEREST!
	!!!==============================================================!!!
	!***CALCULATING MACROSCOPIC VONMISES STRESS***!
	allocate(vm_s(6))	
	s = 1	
	do while(s<=6)
		vm_s(s)=0.0
		s = s+1
	end do
	s = 1
	vm_v = 0.0
	do while (s<=size(atomid_arr))
		vm_s(:)=vm_s(:)+iprop_arr(7:12,s)*iprop_arr(13,s)
		vm_v = vm_v+iprop_arr(13,s)
		s = s+1
	end do
	vm_s(:) =vm_s(:)/vm_v 
	Macro_VM = sqrt((0.5)*((vm_s(1)-vm_s(2))**2+(vm_s(2)-vm_s(3))**2+(vm_s(3)-vm_s(1))**2+6*((vm_s(4))**2+(vm_s(5))**2+(vm_s(6))**2)))
	!!!==============================================================!!!
	!!===================================!!
	!CALCULATING THE AVERAGE PROPERTIES.
	grnmax = maxval(iprop_arr(6,:))
	structmax = maxval(iprop_arr(5,:))
	allocate(avg_props(structmax+4,12,grnmax))	!A 3D BLOCK OF AVERAGE PROPERTIES**12!
	!INITIALISE THE 3D ARRAY TO '0' TO AVOID GARBAGE VALUES
	do n=1,structmax+4		!STRUCTMAX +1 GIVES ALL THE STRUCTURES PRESENT IN THE FILE.  
		do r=1,12	
			do p=1,grnmax
				avg_props(n,r,p)=0
			end do
		end do
	end do
	!!**LEVEL 1 = STRUCTURE TYPE; COLUMNS = AVERAGE PROPERTIES, NUMBER OF ATOMS, GRAIN VOLUME, VON MISES STRESS; ROWS = GRAIN NUMBER**!!
	q = 1
	n = 1
	do q = 1,size(atomid_arr)		
		avg_props(int(iprop_arr(5,q))+1,4:9,int(iprop_arr(6,q)))= &
		& avg_props(int(iprop_arr(5,q))+1,4:9,int(iprop_arr(6,q))) + iprop_arr(7:12,q)*iprop_arr(13,q)	!STRESS * VOLUME
		avg_props(int(iprop_arr(5,q))+1,10,int(iprop_arr(6,q))) = &
		& avg_props(int(iprop_arr(5,q))+1,10,int(iprop_arr(6,q))) + iprop_arr(13,q)		!SUMMATION OF VOLUME FOR A GRAIN
		avg_props(int(iprop_arr(5,q))+1,1:3,int(iprop_arr(6,q))) = &
		& avg_props(int(iprop_arr(5,q))+1,1:3,int(iprop_arr(6,q))) + iprop_arr(2:4,q)	!SUMMATION OF COORDINATES FOR A GRAIN
		avg_props(int(iprop_arr(5,q))+1,11,int(iprop_arr(6,q))) = &
		& avg_props(int(iprop_arr(5,q))+1,11,int(iprop_arr(6,q))) + 1	!NUMBER OF ATOMS IN A PARTICULAR GRAIN FOR A STRUCTURE TYPE 
	end do
	!!**ADDING THREE ELEMENTS IN LEVEL-1: 4->FCC+HCP, 5->FCC+HCP+OTHER DEFECTS, 6->ALL STRUCTURE TYPES**!!
	p = 1
	r = 1
	do p = 1,grnmax
		avg_props(structmax+2,:,p)=avg_props(2,:,p)+avg_props(3,:,p)
		avg_props(structmax+3,:,p)=avg_props(1,:,p)+avg_props(2,:,p)+avg_props(3,:,p)
		do r=1,structmax+1
			avg_props(structmax+4,:,p)=avg_props(structmax+4,:,p)+avg_props(r,:,p)
		end do
	end do
	
	!*3D ARRAY IS OBTAINED CONTAINING ALL THE NECESSARY AVERAGE PROPERTIES*!
	!*VON MISES AND AVERAGE STRESS FOR EACH GRAIN BASED ON STRUCTURE TYPE IS OBTAINED BELOW*!
	p = 1
	r = 1
	do r = 1,structmax+4	
		do p = 1,grnmax
			avg_props(r,4:9,p) = avg_props(r,4:9,p)/avg_props(r,10,p)  !AVERAGE STRESS
			avg_props(r,1:3,p) = avg_props(r,1:3,p)/avg_props(r,11,p)	! GRAIN CENTRE COORDINATES
			!**Von Mises Stress**!
			x = avg_props(r,4,p)-avg_props(r,5,p)
			y = avg_props(r,5,p)-avg_props(r,6,p)
			z = avg_props(r,6,p)-avg_props(r,4,p)
			xy = avg_props(r,9,p)
			yz = avg_props(r,7,p)
			zx = avg_props(r,8,p)
			avg_props(r,12,p) = SQRT((0.5)*(x**2 + y**2 + z**2 )+ 3 *(xy**2+yz**2+zx**2))	
		end do
	end do
	!!!=============================================================================!!!
	!!!===========================================================================!!!
	!!WRITING THE DATA TO A TEXT FILE!!
	!For other defect structure types!
	print*,"DEAR USER, THE AVERAGE PROPERTIES ARE CALCULATED AND BEING WRITTEN INTO A FILE"
	!!!***OPEN THE FILE INTO WHICH THE RESULT IS TO BE STORED***!!!
	open(210, file = 'DATAI_AVGPROP_TESTCASE_4.txt', status='new')  !!!***THE FILE NAME INTO WHICH THE DATA IS TO BE WRITTEN IS SPECIFIED***!!!
	write(210,*) "#*********************PERSONAL PROGRAMMING PROJECT**************************"
	write(210,*) "#AUTHOR: PADMANABHA PAVAN CHANDRA VUNDURTHY(62750) (MSc.Computational Materials Science) &
	& under the supervision of DR. ARUN PRAKASH"
	write(210,*) "#TITLE: Big Data Visualization of Large and Ultra Large Scale&
	& Atomistic Simulation"
	write(210,*) "#Analysis of Average Properties of the Grain separated &
	& in blocks of structure types: 0=Other Defects,1=FCC,2=BCC,3=HCP, 4=BCC..&
	& FCC+HCP, FCC+HCP+0, all structure types."
	do t = 1,num_files
		write(210,*) "#file name->",t,":=",trim(filename_arr(t))
		write(210,*) "#file attribites in sequence::==", file_attributes_arr(:,t)
	end do
	write(210,*) "Macroscopic Von Mises Stress:->", Macro_VM
	u = 1
	o = 1 
	do o=1,grnmax
		summer = summer+avg_props(6,11,o)
	end do
	write(210,*)"Total Number of grains",summer
!	do u = 1,structmax+4
!		do o = 1,grnmax
!			if (o==1)then
!				if(u==structmax+2)then
!					write(210,*) "=======================STRUCTURE TYPES:FCC+HCP======================="
!				else if(u==structmax+3)then
!					write(210,*) "=======================STRUCTURE TYPES:FCC+HCP+other defects======================="
!				else if (u==structmax+4)then
!					write(210,*) "=======================ALL STRUCTURE TYPES======================="
!				else
!					write(210,*) "=======================STRUCTURE TYPE->",u-1,"======================="
!				end if				
!				write(210,*) "GRAIN NUMBER| COORDINATES: X Y Z | STRESS XX YY ZZ ZY XZ XY | GRAIN VOLUME | NUMBER OF ATOMS | VON-MISES STRESS |"
!			end if
!			write(210,*) o,avg_props(u,:,o)
!		end do
!	end do
!	close(210)
	
	print*, "CHECK OUT THE AVERAGE PROPERTIES BY OPENING DATAI_AVGPROP_TESTCASE_4.txt FILE"
	
	!!!**THE ANALYSIS IS COMPLETE*THE RESULTS ARE STORED IN A FILE FOR FURTHER VISUAL ANALYSIS**!!!
	!***END PROGRAM***!
end program Datai_avg
