/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel
 * @Anna 
 * Volume object: This class contains the object and assumes that the distance between the voxels in x,y and z are 1 
 */
public class Volume {
    
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
    public short getVoxel(int x, int y, int z) {
        return data[x + dimX*(y + dimY * z)];
    }
    
    public short getVoxel(int i) {
        return data[i];
    }
    
    public int getDimX() {
        return dimX;
    }
    
    public int getDimY() {
        return dimY;
    }
    
    public int getDimZ() {
        return dimZ;
    }
    
    public void setVoxel(int x, int y, int z, short value) {
        data[x + dimX*(y + dimY*z)] = value;
    }

    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    // -------------------------------------------------------------------------
    
    public short getVoxelNearest(double[] coord) {
        
        // Check for coords outside dimensions
        if (coord[0] < 0 || coord[0] > (dimX-1) || 
            coord[1] < 0 || coord[1] > (dimY-1) || 
            coord[2] < 0 || coord[2] > (dimZ-1) ){
            return 0;
        }
        
        // Find nearest voxel
        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        
        return getVoxel(x, y, z);
    }
    
    /*
    to be implemented: get the trilinear interpolated value. 
    The current implementation gets the Nearest Neightbour
    */
    public short getVoxelInterpolate(double[] coord) {
        
        // Check for coords outside dimensions
        if (coord[0] < 0 || coord[0] > (dimX-1) || 
            coord[1] < 0 || coord[1] > (dimY-1) || 
            coord[2] < 0 || coord[2] > (dimZ-1) ){
            return 0;
        }
        
        // Get i and t for all axis
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        
        double dx = coord[0] - x0;
        double dy = coord[1] - y0;
        double dz = coord[2] - z0;
        
        // Check if i + 1 outside dimensions
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;
        
        if (x1 > dimX-1) x1 = x0;
        if (y1 > dimY-1) y1 = y0;
        if (z1 > dimZ-1) z1 = z0;
        
        // Get surrounding voxels
        short voxel000 = getVoxel(x0, y0, z0);
        short voxel100 = getVoxel(x1, y0, z0);
        short voxel010 = getVoxel(x0, y1, z0);
        short voxel001 = getVoxel(x0, y0, z1);
        short voxel110 = getVoxel(x1, y1, z0);
        short voxel101 = getVoxel(x1, y0, z1);
        short voxel011 = getVoxel(x0, y1, z1);
        short voxel111 = getVoxel(x1, y1, z1);
        
        // Calculate averages of x axis
        short voxelx00 = interpolate(voxel000, voxel100, dx);
        short voxelx01 = interpolate(voxel001, voxel101, dx);
        short voxelx10 = interpolate(voxel010, voxel110, dx);
        short voxelx11 = interpolate(voxel011, voxel111, dx);
        
        // Calculate averages of y axis
        short voxelxy0 = interpolate(voxelx00, voxelx10, dy);
        short voxelxy1 = interpolate(voxelx01, voxelx11, dy);
        
        // Calculate averages of z axis
        short voxelxyz = interpolate(voxelxy0, voxelxy1, dz);
        
        // Return interpolated voxel
        return voxelxyz;
    }
    
    private short interpolate(short v0, short v1, double factor) {
        return (short) Math.round((1.0-factor)*v0 + factor*v1);
    }
    
    // -------------------------------------------------------------------------

    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    public int[] getHistogram() {
        return histogram;
    }
    
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
