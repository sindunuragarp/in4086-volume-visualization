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
    
    public short getVoxelInterpolate(double[] coord) {
        
        // Check for coords outside dimensions
        if (coord[0] < 0 || coord[0] > (dimX-1) || 
            coord[1] < 0 || coord[1] > (dimY-1) || 
            coord[2] < 0 || coord[2] > (dimZ-1) ){
            return 0;
        }
        
        // Get i and t for all axis
        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        double dx = coord[0] - x;
        double dy = coord[1] - y;
        double dz = coord[2] - z;
        
        // Check if i + 1 outside dimensions
        int xi = x + 1;
        int yi = y + 1;
        int zi = z + 1;
        
        if (xi > dimX-1) xi = x;
        if (yi > dimY-1) yi = y;
        if (zi > dimZ-1) zi = z;
        
        // Get surrounding voxels
        short voxel000 = getVoxel(x, y, z);
        short voxel100 = getVoxel(xi, y, z);
        short voxel010 = getVoxel(x, yi, z);
        short voxel001 = getVoxel(x, y, zi);
        short voxel110 = getVoxel(xi, yi, z);
        short voxel101 = getVoxel(xi, y, zi);
        short voxel011 = getVoxel(x, yi, zi);
        short voxel111 = getVoxel(xi, yi, zi);
        
        // Calculate averages of x axis
        short voxelx00 = calcAverageVoxel(voxel000, voxel100, dx);
        short voxelx01 = calcAverageVoxel(voxel001, voxel101, dx);
        short voxelx10 = calcAverageVoxel(voxel010, voxel110, dx);
        short voxelx11 = calcAverageVoxel(voxel011, voxel111, dx);
        
        // Calculate averages of y axis
        short voxelxy0 = calcAverageVoxel(voxelx00, voxelx10, dy);
        short voxelxy1 = calcAverageVoxel(voxelx01, voxelx11, dy);
        
        // Calculate averages of z axis
        short voxelxyz = calcAverageVoxel(voxelxy0, voxelxy1, dz);
        
        // Return interpolated voxel
        return voxelxyz;
    }
    
    private short calcAverageVoxel(short v0, short v1, double dv) {
        return (short) Math.round((1.0-dv)*v0 + dv*v1);
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
