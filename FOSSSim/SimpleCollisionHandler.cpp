#include "SimpleCollisionHandler.h"
#include <stdio.h>
#include <iostream>
#include <set>
#include "DebugHelpers.h"

#define INFINITY std::numeric_limits<double>::infinity()

// BEGIN STUDENT CODE //


// Detects whether two particles are overlapping (including the radii of each)
// and approaching.
// If the two particles overlap and are approaching, returns true and sets 
// the vector n to be the vector between the first and second particle.
// Inputs:
//   scene: The scene data structure. The positions and radii of the particles
//          can be obtained from here.
//   idx1:  The index of the first particle. (Ie, the degrees of freedom
//          corresponding to this particle are entries 2*idx1 and 2*idx1+1 in
//          scene.getX().
//   idx2:  The index of the second particle.
// Outputs:
//   n: The vector between the two particles.
//   Returns true if the two particles overlap and are approaching.
bool SimpleCollisionHandler::detectParticleParticle(TwoDScene &scene, int idx1, int idx2, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*idx1);
    VectorXs x2 = scene.getX().segment<2>(2*idx2);
    
    VectorXs v1 = scene.getV().segment<2>(2*idx1);
    VectorXs v2 = scene.getV().segment<2>(2*idx2);
    
    scalar r1 = scene.getRadius(idx1);
    scalar r2 = scene.getRadius(idx2);
    
    VectorXs movement = x2 - x1;    
    
    VectorXs relativeV = v1 - v2;
    
    n = movement;
    
    #if DEBUG_MODE
        //DEBUGPrintVector(movement);
        //DEBUGPrintVector(n);
    #endif
        
    scalar distance = movement.norm();
    scalar minDistance = r1 + r2;
    
    //movement.normalize();
    //relativeV.normalize();
    
    scalar direction = movement.dot(relativeV);
       
    bool isColliding = false;
    if(direction > 0) {
        isColliding =  distance < minDistance;
    }
    #if DEBUG_MODE
        printf("direction: %.4f, distance: %.4f, min distance: %.4f, isColliding: %d\n", direction, distance, minDistance, isColliding);
    #endif        
    return isColliding;                 
}

// Detects whether a particle and an edge are overlapping (including the radii 
// of both) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the edge.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge. (Ie, the indices of particle with index e are
//          scene.getEdges()[e].first and scene.getEdges()[e].second.)
// Outputs:
//   n: The shortest vector between the particle and the edge.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleEdge(TwoDScene &scene, int vidx, int eidx, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    
    int eVIdx1 = scene.getEdges()[eidx].first;
    int eVIdx2 = scene.getEdges()[eidx].second;
    
    VectorXs x2 = scene.getX().segment<2>(2*eVIdx1);
    VectorXs x3 = scene.getX().segment<2>(2*eVIdx2);      
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*eVIdx1);
    VectorXs v3 = scene.getV().segment<2>(2*eVIdx2);
    
    scalar alpha = (x1 - x2).dot(x3-x2) / ((x3-x2).norm() * (x3-x2).norm());       
    
    if(alpha > 1) {
        alpha = 1;
    } else if(alpha < 0) {
        alpha = 0;
    }
    
    scalar particleRadius = scene.getRadius(vidx);    
    scalar x2Radius = scene.getRadius(scene.getEdges()[eidx].first);
    scalar x3Radius = scene.getRadius(scene.getEdges()[eidx].second);
    
    //scalar edgeRadius = ((x3Radius - x2Radius) * alpha) + x2Radius;
    scalar edgeRadius = scene.getEdgeRadii()[eidx];
    
    VectorXs xAlpha = x2 + (alpha * (x3 - x2));
    VectorXs vAlpha = v2 + (alpha * (v3 - v2));
    
    n = xAlpha - x1;
    
    scalar distance = n.norm();    
    scalar direction = (v1 - vAlpha).dot(n);
                                              
    bool isColliding = false;
    if(direction > 0) {
        if(distance < edgeRadius + particleRadius) {
            isColliding = true;
        }
    }
    
    // your implementation here
    #if DEBUG_MODE
        printf("alpha: %.4f\n, particleRadius: %.4f, edgeRadius: %.4f, direction: %.4f, isColliding: %d\n", alpha, particleRadius, edgeRadius, direction, isColliding);                                 
        //DEBUGPrintVector(n);
    #endif          
    
    return isColliding;
}

// Detects whether a particle and a half-plane are overlapping (including the 
// radius of the particle) and are approaching.
// If the two objects overlap and are approaching, returns true and sets the 
// vector n to be the shortest vector between the particle and the half-plane.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the halfplane. The vectors (px, py) and (nx, ny) can
//          be retrieved by calling scene.getHalfplane(pidx).
// Outputs:
//   n: The shortest vector between the particle and the half-plane.
//   Returns true if the two objects overlap and are approaching.
bool SimpleCollisionHandler::detectParticleHalfplane(TwoDScene &scene, int vidx, int pidx, Vector2s &n)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs px = scene.getHalfplane(pidx).first;
    VectorXs pn = scene.getHalfplane(pidx).second;
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);      
    scalar r1 = scene.getRadius(vidx);
    
    // your implementation here
    n = ((px - x1).dot(pn)/(pn.norm() * pn.norm())) * pn;
    
    scalar distance = n.norm();
    
    bool isColliding = false;    
    if(v1.dot(n) > 0) {
        isColliding = distance < r1;
    }
    
    #if DEBUG_MODE            
        printf("half plane isColliding: %d, n:\n", isColliding);
        DEBUGPrintVector(n);
    #endif
    
    return isColliding;
}


// Responds to a collision detected between two particles by applying an impulse
// to the velocities of each one.
// You can get the COR of the simulation by calling getCOR().
// Inputs:
//   scene: The scene data structure.
//   idx1:  The index of the first particle.
//   idx2:  The index of the second particle.
//   n:     The vector between the first and second particle.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleParticle(TwoDScene &scene, int idx1, int idx2, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    VectorXs &v = scene.getV();
    
    // your implementation here    
    // 
    VectorXs nhat = n;
    nhat.normalize();
    
    VectorXs v1 = v.segment<2>(2*idx1);
    VectorXs v2 = v.segment<2>(2*idx2);    
    
    bool pIsFixed1 = scene.isFixed(idx1);
    bool pIsFixed2 = scene.isFixed(idx2);
    
    scalar m1 = pIsFixed1 ? INFINITY : M[2*idx1];
    scalar m2 = pIsFixed2 ? INFINITY : M[2*idx2];
    
    scalar corFactor = (1 + getCOR()) / 2;        
    
    #if DEBUG_MODE
        //printf("----------------------\n");
        //printf("mass vectors\n");
        //DEBUGPrintVector(m1);
        //DEBUGPrintVector(m2);
        printf("cor: %.4f, corfactor: %.4f", getCOR(), corFactor);
    #endif     
    
    if(!pIsFixed1) {
        VectorXs i1 = corFactor * (((2*(v2-v1)).dot(nhat)) / (1 + (m1/m2))) * nhat;    
        v.segment<2>(2*idx1) += i1;
    }
    
    if(!pIsFixed2) {
        VectorXs i2 = corFactor * (((2*(v2-v1)).dot(nhat)) / ((m2/m1) + 1)) * nhat;    
        v.segment<2>(2*idx2) -= i2;
    }               
}

// Responds to a collision detected between a particle and an edge by applying
// an impulse to the velocities of each one.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   eidx:  The index of the edge.
//   n:     The shortest vector between the particle and the edge.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleEdge(TwoDScene &scene, int vidx, int eidx, const Vector2s &n)
{
    const VectorXs &M = scene.getM();
    
    VectorXs nhat = n;
    nhat.normalize();
    
    int eidx1 = scene.getEdges()[eidx].first;
    int eidx2 = scene.getEdges()[eidx].second;
    
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*eidx1);
    VectorXs x3 = scene.getX().segment<2>(2*eidx2);
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*eidx1);
    VectorXs v3 = scene.getV().segment<2>(2*eidx2);
    
    bool pIsFixed1 = scene.isFixed(vidx);
    bool pIsFixed2 = scene.isFixed(eidx1);
    bool pIsFixed3 = scene.isFixed(eidx2);
    
    scalar m1 = pIsFixed1 ? INFINITY : M[2*vidx];
    scalar m2 = pIsFixed2 ? INFINITY : M[2*eidx1];
    scalar m3 = pIsFixed3 ? INFINITY: M[2*eidx2];
    
    // your implementation here
    // 
    scalar alpha = (x1 - x2).dot(x3-x2) / ((x3-x2).norm() * (x3-x2).norm());       
    
    if(alpha > 1) {
        alpha = 1;
    } else if(alpha < 0) {
        alpha = 0;
    }
    
    scalar restAlpha = 1 - alpha;
    
    VectorXs ve = (restAlpha * v2) + ((alpha) * v3);   
    
    scalar corFactor = (1 + getCOR()) / 2;
    
    #if DEBUG_MODE 
        printf("cor: %.4f, corfactor: %.4f", getCOR(), corFactor);
    #endif
    
    if(!pIsFixed1) {
        VectorXs i1 = corFactor * (((2 * (ve - v1)).dot(nhat)) / (1 + ((restAlpha * restAlpha * m1)/m2) + ((alpha * alpha * m1) / m3))) * nhat;        
        scene.getV().segment<2>(2*vidx) += i1;
    }
    
    if(!pIsFixed2) {
        VectorXs i2 = corFactor * (((2 * restAlpha * (ve - v1)).dot(nhat))/((m2/m1) + (restAlpha * restAlpha) + ((alpha * alpha * m2)/m3))) * nhat;    
        scene.getV().segment<2>(2*eidx1) -= i2;        
    }
    
    if(!pIsFixed3) {
        VectorXs i3 = corFactor * ((((2 * alpha)*(ve - v1)).dot(nhat)) / ((m3/m1) + ((restAlpha * restAlpha * m3)/m2) + (alpha * alpha))) * nhat;        
        scene.getV().segment<2>(2*eidx2) -= i3;
    }
        
}


// Responds to a collision detected between a particle and a half-plane by 
// applying an impulse to the velocity of the particle.
// Inputs:
//   scene: The scene data structure.
//   vidx:  The index of the particle.
//   pidx:  The index of the half-plane.
//   n:     The shortest vector between the particle and the half-plane.
// Outputs:
//   None.
void SimpleCollisionHandler::respondParticleHalfplane(TwoDScene &scene, int vidx, int pidx, const Vector2s &n)
{
    VectorXs nhat = n;
    
    nhat.normalize();    
    
    VectorXs v = scene.getV().segment<2>(2*vidx);
   
    scalar corFactor = (1.0 + getCOR()) / 2.0;   
    
    VectorXs i = corFactor * (2.0 * ((v.dot(nhat)) * nhat));
    
    #if DEBUG_MODE 
        printf("cor: %.4f, corfactor: %.4f, impulse: \n", getCOR(), corFactor);
        DEBUGPrintVector(i);
        printf("n:\n");
        DEBUGPrintVector(n);
        printf("nhat:\n");
        DEBUGPrintVector(nhat);
    #endif      
      
    scene.getV().segment<2>(2*vidx) -= i;
    
    // your implementation here
    
}
