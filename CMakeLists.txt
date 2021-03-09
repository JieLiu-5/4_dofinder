cmake_minimum_required(VERSION 3.10.2)
project (dofinder)                                           

#Bring the headers, such as Student.h into the project
include_directories(include)


#Can manually add the sources using the set command as follows:
# set(SOURCES src/test.cpp src/functions.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")                                     # we look for .cpp file automatically based on the main() function.
                                                                   # therefore, only one main() function is allowed

                                                                   

message( "OOOOOOOOOOOOO OOOOOOOOOOOOOO OOOOOOOOOOOOOOOOOOOOOOOOOO OOOOOOOOOOOOOO OOOOOOOOOOOOO" )
message( "OOOOOOOOOOOO   OOOOOOOOOOOO   OOOOOOOOOOOOOOOOOOOOOOOO   OOOOOOOOOOOO   OOOOOOOOOOOO" )
message( "OOOOOOOOOOO     OOOOOOOOOO     OOOOOOOOOOOOOOOOOOOOOO     OOOOOOOOOO     OOOOOOOOOOO" )
message( "OOOOOOOOOO   O   OOOOOOOO   O   OOOOOOOOOOOOOOOOOOOO   O   OOOOOOOO   O   OOOOOOOOOO" )
message( "OOOOOOOOO   O O   OOOOOO   O O   OOOOOOOOOOOOOOOOOO   O O   OOOOOO   O O   OOOOOOOOO" )
message( "OOOOOOOO   O   O   OOOO   O   O   OOOOOOOOOOOOOOOO   O   O   OOOO   O   O   OOOOOOOO" )
message( "OOOOOOOO   O | O   OOOO   O | O   OOOOOOOOOOOOOOOO   O | O   OOOO   O | O   OOOOOOOO" )
message( "OOOOOOOOO    O    OOOOOO    O    OOOOOOOOOOOOOOOOOO    O    OOOOOO    O    OOOOOOOOO" )
message( "OOOOOOOOOO       OOOOOOOO       OOOOOOOOOOOOOOOOOOOO       OOOOOOOO       OOOOOOOOOO" )
message( "OOOOOOOOOOOO   OOOOOOOOOOOO   OOOOOOOOOOOOOOOOOOOOOOOO   OOOOOOOOOOOO   OOOOOOOOOOOO" )
message( "OOOOOOOOOOOOO OOOOOOOOOOOOOO OOOOOOOOOOOOOOOOOOOOOOOOOO OOOOOOOOOOOOOO OOOOOOOOOOOOO" )
message( "OOOOOOOOOOOOO OOOOOOOOOOOOOO OOOOOOOOOOOOOOOOOOOOOOOOOO OOOOOOOOOOOOOO OOOOOOOOOOOOO" )                                                             

add_executable(dofinder ${SOURCES})                             # this determines the name of the generated executable program
