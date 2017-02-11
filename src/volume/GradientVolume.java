/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 * @ Anna
 * This class contains the pre-computes gradients of the volume. This means calculates the gradient
 * at all voxel positions, and provides functions
 * to get the gradient at any position in the volume also continuous..
*/
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        grad1 = new VoxelGradient[dimX * dimY * dimZ];
        grad2 = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
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
    
    public VoxelGradient getVoxel(int i) {
        return grad1[i];
    }
    
    public VoxelGradient getGradient(int x, int y, int z) {
        return getGradient(x,y,z,false);
    }
    
    public VoxelGradient getGradient(int x, int y, int z, boolean isSecondGradient) {
        if(isSecondGradient) return grad2[x + dimX * (y + dimY * z)];
        return grad1[x + dimX * (y + dimY * z)];
    }
    
    public VoxelGradient getGradient(double[] coord) {
        return getGradientInterpolate(coord);
    }
    
    public VoxelGradient getSecondGradient(double[] coord) {
        return getGradientInterpolate(coord, true);
    }
    
    public void setVoxel(int i, VoxelGradient value) {
        grad1[i] = value;
    }
    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        setGradient(x, y, z, value, false);
    }
    
    public void setGradient(int x, int y, int z, VoxelGradient value, boolean isSecondGradient) {
        if(isSecondGradient) grad2[x + dimX * (y + dimY * z)] = value;
        else grad1[x + dimX * (y + dimY * z)] = value;
    }

    // -------------------------------------------------------------------------
    
    public VoxelGradient getGradientNearest(double[] coord) {
        return getGradientNearest(coord, false);
    }
    
    public VoxelGradient getGradientNearest(double[] coord, boolean isSecondGradient) {
        
        // Check for coords outside dimensions
        if (coord[0] < 0 || coord[0] > (dimX-1) || 
            coord[1] < 0 || coord[1] > (dimY-1) || 
            coord[2] < 0 || coord[2] > (dimZ-1) ){
            return zero;
        }

        // Find nearest voxel
        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        
        return getGradient(x, y, z, isSecondGradient);
    }

    /*
    To be implemented: Returns trilinear interpolated gradient based on the precomputed gradients. 
    Use function interpolate. Use getGradientNN as bases
    */
    public VoxelGradient getGradientInterpolate(double[] coord) {
        return getGradientInterpolate(coord, false);
    }
    
    public VoxelGradient getGradientInterpolate(double[] coord, boolean isSecondGradient) {

        // Check for coords outside dimensions
        if (coord[0] < 0 || coord[0] > (dimX-1) || 
            coord[1] < 0 || coord[1] > (dimY-1) || 
            coord[2] < 0 || coord[2] > (dimZ-1) ){
            return zero;
        }
        
        // Get i and t for all axis
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        
        float dx = (float)(coord[0] - x0);
        float dy = (float)(coord[1] - y0);
        float dz = (float)(coord[2] - z0);
        
        // Check if i + 1 outside dimensions
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;
        
        if (x1 > dimX-1) x1 = x0;
        if (y1 > dimY-1) y1 = y0;
        if (z1 > dimZ-1) z1 = z0;
        
        // Get surrounding grads
        VoxelGradient grad000 = getGradient(x0, y0, z0, isSecondGradient);
        VoxelGradient grad100 = getGradient(x1, y0, z0, isSecondGradient);
        VoxelGradient grad010 = getGradient(x0, y1, z0, isSecondGradient);
        VoxelGradient grad001 = getGradient(x0, y0, z1, isSecondGradient);
        VoxelGradient grad110 = getGradient(x1, y1, z0, isSecondGradient);
        VoxelGradient grad101 = getGradient(x1, y0, z1, isSecondGradient);
        VoxelGradient grad011 = getGradient(x0, y1, z1, isSecondGradient);
        VoxelGradient grad111 = getGradient(x1, y1, z1, isSecondGradient);
        
        // Calculate averages of x axis
        VoxelGradient gradx00 = interpolate(grad000, grad100, dx);
        VoxelGradient gradx01 = interpolate(grad001, grad101, dx);
        VoxelGradient gradx10 = interpolate(grad010, grad110, dx);
        VoxelGradient gradx11 = interpolate(grad011, grad111, dx);
        
        // Calculate averages of y axis
        VoxelGradient gradxy0 = interpolate(gradx00, gradx10, dy);
        VoxelGradient gradxy1 = interpolate(gradx01, gradx11, dy);
        
        // Calculate averages of z axis
        VoxelGradient gradxyz = interpolate(gradxy0, gradxy1, dz);
        
        // Return interpolated grad
        return gradxyz;
    }
    
    /*
    To be implemented: this function linearly interpolates gradient vector g0 and g1 given the factor (t) 
    the resut is given at result. You can use it to tri-linearly interpolate the gradient  
    */
    private VoxelGradient interpolate(VoxelGradient g0, VoxelGradient g1, float factor) {
        
        float x = Math.round((1.0-factor)*g0.x + factor*g1.x);
        float y = Math.round((1.0-factor)*g0.y + factor*g1.y);
        float z = Math.round((1.0-factor)*g0.z + factor*g1.z);
        
        return new VoxelGradient(x,y,z);
    }
    
    // -------------------------------------------------------------------------
    
    /*
    To be implemented: compute the gradient of contained in the volume attribute
    */
    private void compute() {
        // clear data
        for (int i=0; i<grad1.length; i++) {
            grad1[i] = zero;
        }
        
        // calculate gradients
        for(int x=1; x<dimX-1; x++) {
            for(int y=1; y<dimY-1; y++) {
                for(int z=1; z<dimZ-1; z++) {
                    
                    // get voxels
                    short o  = volume.getVoxel(x, y, z);
                    short x0 = volume.getVoxel(x-1, y, z);
                    short x1 = volume.getVoxel(x+1, y, z);
                    short y0 = volume.getVoxel(x, y-1, z);
                    short y1 = volume.getVoxel(x, y+1, z);
                    short z0 = volume.getVoxel(x, y, z-1);
                    short z1 = volume.getVoxel(x, y, z+1);
                    
                    // calculate central differences
                    float dx = (float) (x1-x0)/2;
                    float dy = (float) (y1-y0)/2;
                    float dz = (float) (z1-z0)/2;
                    
                    // set gradient
                    setGradient(x, y, z, new VoxelGradient(dx, dy, dz), false);
                    
                    // calculate second order central differences
                    float d2x = (float) (x1+2*o-x0)/4;
                    float d2y = (float) (y1+2*o-y0)/4;
                    float d2z = (float) (z1+2*o-z0)/4;
                    
                    // set second gradient
                    setGradient(x, y, z, new VoxelGradient(d2x, d2y, d2z), true);
                }
            }
        }
    }
    
    /* 
    to be implemented: Returns the maximum gradient magnitude
    */
    public float getMaxGradientMagnitude() {
        float maximum = grad1[0].mag;
        for (VoxelGradient grad : grad1) {
            float now = grad.mag;
            if(now > maximum) maximum = now;
        }
        return maximum;
    }
    
    // -------------------------------------------------------------------------
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] grad1;
    VoxelGradient[] grad2;
    Volume volume;
    double maxmag;
}
