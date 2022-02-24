classdef sC
% Class that serves as "library" (pseudo screwCalculus) for different
% functions useful to compute the kinematics of a shaping machine for
% spiral bevel wheels.
%
%       FUNCTIONS (callable as METHODS: sC.method() )
%
%   N.B. You can improve calling speed of library functions renaming the
%   class in the workspace (i.e. SC = sC; SC.hat() calls the hat() function)
%
% - hat(vector)
%
% - rotX(alpha) (the MATLAB one does not allow for symbolic computations)
%
% - rotY(alpha)
%
% - rotZ(alpha)
%
% - hom_mat(rotationMatrix, displacement)
%
% - hom_mat_inv(Gst): inverse of a homog. matrix transf. Gst
%
% - adjoint(rotationMatrix,displacement) or adjoint(homogeneousmatrix)
%
% - expSkew(axis,th) -> exp(omegahat*theta)
%
% - unitTwist(jointType,axis,q) (NO HELIC. JOINTS): reference unit Twist
%
% - expTw(unitTwist,th) ( NO HELIC. JOINS): twist Exponential
%
% - FWKin(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn}): forward Kinematics
%
% - AD_FWKin(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn}): automatic diff.
%   and Forward Kinematics (automatic diff. variable on joint vars).
%
% - SpJac({Y1,var1},{Y2,var2},...,{Yn,varn}): spatial Jacobian
%
% - BodyJac(gst0,{Y1,var1},{Y2,var2},...,{Yn,varn}): body Jacobian
%
%   TODO: copy external help function or what?

    methods(Static)
        % Function signatures only.
        % Definitions are in a different file for each function
        
        X = hat(v)
        
        [n, theta] = rotToAxisAngle(R)
        
        R = rotX(alpha)
        
        R = rotY(alpha)
        
        R = rotZ(alpha)
        
        R = rotZ2D(alpha)
        
        R = diffRotZ2D(alpha)
        
        Gst = hom_mat(rotm,disp)

        GstInv = hom_mat_inv(Gst)

        GstInv = hom_mat_inv2D(Gst)
        
        X = adjoint(rotm,disp)
        
        X = expSkew(axis,th)
        
        uT = unitTwist(jointType,axis,q)
        
        E = expTw(unitTwist,th)
        
        gst = FWKin(gst0,varargin)
        
        gst = AD_FWKin(gst0,varargin)
        
        Js = SpJac(varargin)
        
        Jb = BodyJac(gst0,varargin)
        
    end
end