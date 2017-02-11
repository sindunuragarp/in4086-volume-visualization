/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;
import java.lang.Thread;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author michel
 * @Anna
 * This class has the main code that generates the raycasting result image. 
 * The connection with the interface is already given.  
 * The different modes mipMode, slicerMode, etc. are already correctly updated
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }
    
    // -------------------------------------------------------------------------

    /* 
    To be partially implemented: This function traces the rays through the 
    volume. Have a look and check that you understand how it works.
    You need to introduce here the different modalities MIP/Compositing/TF2/ etc...
    */
    void raycast(double[] viewMatrix) {
        
        // rendering vars
        int increment = 1;
        float sampleStep = 0.2f;
        
        if(interactiveMode) {
            increment = (int) Math.floor(image.getWidth()/100);
            sampleStep = 1.0f;
        }
        
        // vector uVec and vVec define a plane through the origin, perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // init vars
        int imageCenter = image.getWidth() / 2;
        short maxIntensity = volume.getMaximum();
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        
        // clear image
        for (int i = 0; i < image.getWidth(); i++) {
            for (int j = 0; j < image.getHeight(); j++) {
                image.setRGB(i, j, 0);
            }
        }

        // raycast loop
        ArrayList<Thread> threads = new ArrayList();
        for (int j = 0; j < image.getHeight(); j += increment) {
            for (int i = 0; i < image.getWidth(); i += increment) {
                
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) - viewVec[0] * imageCenter + volume.getDimX() / 2.0;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) - viewVec[1] * imageCenter + volume.getDimY() / 2.0;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) - viewVec[2] * imageCenter + volume.getDimZ() / 2.0;

                computeEntryAndExit(pixelCoord, viewVec, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    
                    // Start calculation in different threads
                    Thread newThread = new ColorSetter(i, j, viewVec, entryPoint, exitPoint, sampleStep, increment, maxIntensity);
                    newThread.run();
                    threads.add(newThread);
                }
            }
        }
        
        // Wait for all threads to terminate
        for (Thread thread : threads) {
            try { thread.join(); } 
            catch (InterruptedException ex) {}
        }
    }
    
    // -------------------------------------------------------------------------
    
    // set shading parameters
    private double specular = 0.2;
    private double diffuse = 0.7;
    private double ambient = 0.1;
    private double[] lightVec = {1.0, 0.0, 0.0};
            
    private class ColorSetter extends Thread {
        
        private final int posI;
        private final int posJ;
        
        private final double[] viewVec;
        private final double[] entryPoint;
        private final double[] exitPoint;
        
        float sampleStep;
        int increment;
        short maxIntensity;
        
        ColorSetter(int i, int j, double[] viewVec, double[] entryPoint, double[] exitPoint, float sampleStep, int increment, short maxIntensity) {
            this.posI = i;
            this.posJ = j;
            this.viewVec = viewVec;
            this.entryPoint = entryPoint;
            this.exitPoint = exitPoint;
            this.sampleStep = sampleStep;
            this.increment = increment;
            this.maxIntensity = maxIntensity;
        }

        @Override
        public void run() {
            
            int pixelColor = 0;

            if(mipMode) pixelColor = traceRayMIP(entryPoint,exitPoint);
            else if(compositingMode || tf2dMode) pixelColor = traceRayCompositing(entryPoint,exitPoint);

            // Check for out of bounds
            int ix = posI + increment;
            int jx = posJ + increment;
            
            if (increment > 1) {
                if (ix > image.getWidth()) ix = image.getWidth();
                if (jx > image.getHeight()) jx = image.getHeight();
            }
            
            // Set pixels
            for (int ii = posI; ii < ix; ii++) {
                for (int jj = posJ; jj < jx; jj++) {
                    image.setRGB(ii, jj, pixelColor);
                }
            }
        }
        
        // ---------------------------------------------------------------------
        
        /*
        to be implemented:  You need to sample the ray and implement the MIP
        right now it just returns yellow as a color
        */
        private int traceRayMIP(double[] entryPoint, double[] exitPoint) {
            int totalSteps = calcTotalSteps(entryPoint, exitPoint);
            short voxelMax = 0;

            // find maximum intensity
            for (int i=0; i<totalSteps; i++){
                
                // check max
                double[] coord = calcCoord(i, totalSteps);
                short voxelNow = interactiveMode ? volume.getVoxelNearest(coord) : volume.getVoxelInterpolate(coord);
                if (voxelNow > voxelMax) voxelMax = voxelNow;
            }

            // calculate rgb value
            int a = (int)(255 * voxelMax / maxIntensity);
            int r = 255;
            int g = 255;
            int b = 0;
            
            int argb = (a << 24) | (r << 16) | (g << 8) | b;
            return argb;
        }

        private int traceRayCompositing(double[] entryPoint, double[] exitPoint) {
            int totalSteps = calcTotalSteps(entryPoint, exitPoint);
            TFColor color = new TFColor(0,0,0,0);
            
            // run compositing (front to back)
            for (int i=0; i<totalSteps; i++){
                
                // get coord values
                double[] coord = calcCoord(i, totalSteps);
                VoxelGradient grad = interactiveMode ? gradients.getGradientNearest(coord) : gradients.getGradientInterpolate(coord);
                short voxel = interactiveMode ? volume.getVoxelNearest(coord) : volume.getVoxelInterpolate(coord);
                
                // get voxel color
                TFColor vColor = tf2dMode ? getTF2DColor(voxel, grad) : tFunc.getColor(voxel);
                if(shadingMode) vColor = calcShade(vColor, grad);
                
                // accumulate color (compositing)
                color.r += (1-color.a) * vColor.r * vColor.a;
                color.g += (1-color.a) * vColor.g * vColor.a;
                color.b += (1-color.a) * vColor.b * vColor.a;
                color.a += (1-color.a) * vColor.a;
                
                // early termination
                if(color.a > 0.99) break;
            }
            
            // calculate rgb value
            int a = color.a < 1.0 ? (int) Math.floor(color.a * 255) : 255;
            int r = color.r < 1.0 ? (int) Math.floor(color.r * 255) : 255;
            int g = color.g < 1.0 ? (int) Math.floor(color.g * 255) : 255;
            int b = color.b < 1.0 ? (int) Math.floor(color.b * 255) : 255;

            int argb = (a << 24) | (r << 16) | (g << 8) | b;
            return argb;
        }
        
        // ---------------------------------------------------------------------
        
        private int calcTotalSteps(double[] entryPoint, double[] exitPoint) {
            double xDist = exitPoint[0] - entryPoint[0];
            double yDist = exitPoint[1] - entryPoint[1];
            double zDist = exitPoint[2] - entryPoint[2];

            double rayLength = Math.sqrt(xDist*xDist + yDist*yDist + zDist*zDist);
            int totalSteps = (int) Math.floor(rayLength/sampleStep);
            
            return totalSteps;
        }
        
        private double[] calcCoord(double currentStep, double totalStep) {
            double step = sampleStep * currentStep;
            double[] coord = new double[3];
            
            coord[0] = entryPoint[0] - (step * viewVec[0]);
            coord[1] = entryPoint[1] - (step * viewVec[1]);
            coord[2] = entryPoint[2] - (step * viewVec[2]);
            
            return coord;
        }
        
        // ---------------------------------------------------------------------
        
        private TFColor getTF2DColor(short intensity, VoxelGradient gradient) {
            
            // transfer function values
            TransferFunction2DEditor.TriangleWidget widget = tfEditor2D.triangleWidget;
            TFColor color = widget.color;
            short baseIntensity = widget.baseIntensity;
            double radius = widget.radius;
            
            // get gradient based opacity
            double alpha = calcGradientBasedOpacity(intensity, baseIntensity, radius, gradient.mag);

            // return color
            return new TFColor(color.r, color.g, color.b, color.a * alpha);
        }
        
        private double calcGradientBasedOpacity(short in, short bin, double rad, float mag) {
            
            int din = in - bin;
            
            // absolute opacity conditions
            if(mag == 0 && din == 0) return 1.0;
            if(mag < 0 || Math.abs(din) > rad*mag) return 0.0;
            
            // avoid division by zero
            if(mag == 0) mag = 0.00001f;
            if(rad == 0) rad = 0.00001f;
            
            // calculate partial opacity
            double alpha = 1 - (1/rad)*Math.abs(din/mag);
            return alpha;
        }
        
        private TFColor calcShade(TFColor color, VoxelGradient grad) {
            
            // calculate normalized gradient vector
            double[] gradVec = new double[3];
            gradVec[0] = -grad.x / grad.mag;
            gradVec[1] = -grad.y / grad.mag;
            gradVec[2] = -grad.z / grad.mag;
            
            // TODO : calculate lighting effects
            
            // fuse lighting effects
            double contrast = 1.0;
            return new TFColor(color.r*contrast, color.g*contrast, color.b*contrast, color.a);
        }
    }
    
    // -------------------------------------------------------------------------

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                
                // Get coords
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];

                // Get interpolated voxel
                int val = volume.getVoxelInterpolate(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                // Alternatively, apply the transfer function to obtain a color
                /*
                TFColor auxColor = new TFColor();
                auxColor = tFunc.getColor(val);
                voxelColor.r=auxColor.r;
                voxelColor.g=auxColor.g;
                voxelColor.b=auxColor.b;
                voxelColor.a=auxColor.a;
                */
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    // -------------------------------------------------------------------------
    
    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     
    public void setShadingMode(boolean mode) {
        shadingMode = mode;
        changed();
    }
    
    public void setMIPMode() {
        setMode(false, true, false, false);
    }
    
    public void setSlicerMode() {
        setMode(true, false, false, false);
    }
    
    public void setCompositingMode() {
        setMode(false, false, true, false);
    }
    
    public void setTF2DMode() {
        setMode(false, false, false, true);
    }
    
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        changed();
    }
        
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();
    }
    
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }

    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));
    }

    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            //System.out.println("Plane pos: " + plane_pos[0] + " " + plane_pos[1] + " " + plane_pos[2]);
            //System.out.println("Intersection: " + intersection[0] + " " + intersection[1] + " " + intersection[2]);
            //System.out.println("line_dir * intersection: " + VectorMath.dotproduct(line_dir, plane_normal));

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) > 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }

    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);
    }
    
    // -------------------------------------------------------------------------
    
    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
