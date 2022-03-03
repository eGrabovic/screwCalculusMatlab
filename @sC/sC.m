classdef sC
% Class that serves as "library" (pseudo screwCalculus) for different
% functions useful to compute the kinematics of a shaping machine for
% spiral bevel wheels.
%
%       FUNCTIONS (callable as METHODS: sC.method() )
%
%   N.B. You can improve calling speed of library functions renaming the
%   class in the workspace (i.e. SC = sC; SC.hat() calls the hat() function)

    methods(Static)
        % Function signatures only.
        % Definitions are in a different file for each function
        
        gst = AD_FWKin(screwObj, gst0,varargin)
        
        X = adjoint(screwObj, rotm,disp)
        
        Jb = BodyJac(screwObj, gst0, varargin)
        
        X = expSkew(screwObj, axis,th)
        
        E = expTw(screwObj, unitTwist,th)
        
        gst = FWKin(screwObj, gst0,varargin)
        
        X = hat(screwObj, v)
        
        Gst = hom_mat(screwObj, rotm,disp)

        GstInv = hom_mat_inv(screwObj, Gst)

        GstInv = hom_mat_inv2D(screwObj, Gst)
        
        [n, theta] = rotToAxisAngle(screwObj, R)
        
        R = rotX(screwObj, alpha)
        
        R = rotY(screwObj, alpha)
        
        R = rotZ(screwObj, alpha)
        
        R = rotZ2D(screwObj, alpha)
        
        Js = SpJac(screwObj, varargin)
        
        uT = unitTwist(screwObj, jointType, axis, q)
        
    end
end