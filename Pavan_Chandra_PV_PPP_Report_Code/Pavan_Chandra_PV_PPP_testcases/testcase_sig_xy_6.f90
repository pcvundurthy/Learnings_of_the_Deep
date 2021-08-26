program testcase_sig_xy_6
implicit none
	integer::error, num_files,a=0, num_props, b=0, i,j, errro,k,l,m
	character(120), dimension(:), allocatable:: filename_arr, prop_name_arr
	integer,dimension(:,:), allocatable::file_attributes_arr, prop_attributes_arr
	integer::c, atomid_filenum, d=1, e, f, g, h, span, counter, grnmax, structmax
	integer::n,o, p,q,r,s,t,u,x,y,z,xy,yz,zx,vm_v, Macro_VM
	real, dimension(:,:), allocatable::iprop_arr, temp_arr
	real,dimension(:), allocatable::atomid_arr, vm_s
	real, dimension(:,:,:), allocatable::avg_props
	!!================================!!	
	!OPEN AND READ THE PARAMETER FILE FOR USER INPUT
	open (unit=201, file='Parameter_File_BigData.txt', iostat = error)
	if (error .ne. 0) then
            write(*,*) "file cannot be opened"
    end if
    read(201, *)
    read(201,*)num_files
	print*,"=========",num_files
    read(201, *)
    allocate(filename_arr(num_files))
    allocate(file_attributes_arr(3,num_files))
    !READ FOR STORING THE FILE NAMES AND FILE ATTRIBUTES
    do while(a<num_files)
		read(201,*)filename_arr(a+1), file_attributes_arr(:,a+1)  !!THE INDEX 'a' STARTS AT 0
		a = a+1
	end do
	print*,"=========",filename_arr, file_attributes_arr
	read(201,*)
	read(201,*)
	read(201,*)num_props
	allocate(prop_name_arr(num_props))
	allocate(prop_attributes_arr(3,num_props))
	do while(b<num_props)
		read(201,*)prop_name_arr(b+1), prop_attributes_arr(:,b+1) !!THE INDEX 'b' STARTS AT 0
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
			i = maxval(file_attributes_arr(2,:))
			j = maxval(file_attributes_arr(1,:))-minval(file_attributes_arr(3,:))
			allocate(temp_arr(i,j))
			allocate(atomid_arr(file_attributes_arr(1,atomid_filenum)-file_attributes_arr(3,atomid_filenum)))
			allocate(iprop_arr(13,size(atomid_arr)))
			!!==OPEN THE FILE==!!
			open(unit=202,file = filename_arr(atomid_filenum), iostat = errro)	
			if (errro .ne. 0) then
					write(*,*) "file cannot be opened"
			end if	
			!***==Initializing the Array to '0'==***!
			l=1
			do while(l<=j)
				temp_arr(:,l)=0
				l=l+1
			end do
			!----------------------!
			!***reading the file into a temporary array: temp_arr***!
			k = 1
			do
				if (k<=file_attributes_arr(3,atomid_filenum))then
					read (202,*,end=10)
				else if (k<=file_attributes_arr(1,atomid_filenum)) then
					read (202,*,end=10) temp_arr(1:file_attributes_arr(2,atomid_filenum),k-file_attributes_arr(3,atomid_filenum))					
				end if
				k=k+1				
			end do
			!!----------------------!
10 continue
			close(202)
			!!==CLOSE THE FILE==!!
			!!!-----ASSIGNING TO ATOM ID ARRAY AND THEN ASSIGNING ATOM ID ARRAY TO THE FIRST COLUMN OF THE IPROP ARRAY-----!!!
			atomid_arr = temp_arr(1,1:file_attributes_arr(1,atomid_filenum)-file_attributes_arr(3,atomid_filenum))
			iprop_arr(1,:)=atomid_arr
			!!-----------------------!
			do while (d<=num_props)
				if (prop_attributes_arr(2,d)==atomid_filenum)then					
			!**-------------------------------------------------------**!
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
			!**------------------------------------------------------**!
				end if
				d = d+1
			end do
		end if
		c = c+1
	end do
	!CHECK FOR ATOM ID FILE IS ACCOMPLISHED AND PROPERTIES IN THE SAME FILE ARE EXTRACTED AND STORED IN THE ARRAY.
	!NOW, EXTRACT NECESSARY DATA FROM OTHER FILES IGNORING THE ATOM ID FILE
	e = 1
	do while (e<=num_files)
		if (e==atomid_filenum)then
			goto 25
		end if
		!!!!------Initializing the temporary array again back to '0'----!!!!
		l=1
		do while(l<=j)
			temp_arr(:,l)=0
			l=l+1
		end do
		!!--------------!!
		span = file_attributes_arr(1,e)-file_attributes_arr(3,e)
		!!OPEN A NEW FILE!!
		open(unit=203,file = filename_arr(e), iostat = errro)	
		if (errro .ne. 0) then
				write(*,*) "file cannot be opened"
		end if	
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
		!!CLOSE THE FILE!!
		f = 1
		do while(f<=num_props)			
			if (prop_attributes_arr(2,f)==e)then		
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
25 continue
		e = e+1	
	end do
	!!!Creating a text file for test case!!!
	!!!ASSIGN THE STRESS VALUES FOR THE ARRAY!!!
	grnmax = maxval(iprop_arr(6,:))
	structmax = maxval(iprop_arr(5,:))
	do t = 1,int(grnmax)
		do u=1,size(atomid_arr)
			if (iprop_arr(6,u)==t) then
				iprop_arr(7:12,u)=(/ 0.0, 0.0, 0.0, 0.0, 0.0, 6.0*t/)	!!**sigma_xx=1.0*grain num, and other stresses are 0**!!
				iprop_arr(13,u)=0.5	!!**the volume is 0.5 and remains constain for each grain**!!
			end if
		end do
	end do	
	!!** OPEN THE FILE FOR WRITING THE DATA**!!
	open(204, file = 'DATA_FOR_TESTCASE_6.txt', status='new') 
	write(204,*)"#Personal Programming Project---Pavan Chandra Vundurthy(62750)---under the supervision of Dr. Arun Prakash"
	write(204,*)"#Below is the Generated data for the test case_6 with corresponding stresses to check for the analysis"
	write(204,*)"=============================================================================================="""
	o = 1
	do o=1,size(atomid_arr)		
		write(204,*) iprop_arr(:,o)
	end do
	close(204)
	!!**CLOSED THE FILE**!!

end program testcase_sig_xy_6
