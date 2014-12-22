/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package canvas;

import control.Plotter;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.MouseEvent;
import javax.swing.event.MouseInputListener;

/**
 *
 * @author edsonlan
 */
public class NativeCanvas extends Canvas implements MouseInputListener,ComponentListener {

 
    private int id_;
    private Plotter plotter_;

    public NativeCanvas(int id,Plotter plotter) {
        plotter_=plotter;
        id_=id;
        addMouseListener(this);
        addMouseMotionListener(this);
//        addComponentListener(this);
        setBackground(Color.black);

    }

    public native void paint(Graphics g);

    @Override
    public void mouseClicked(MouseEvent e) {

    }

    @Override
    public void mousePressed(MouseEvent e) {
//        plotter_.cleanList();
        plotter_.plot("level", e.getX(),e.getY());
        paint(getGraphics());
        plotter_.update();
    }

    @Override
    public void mouseReleased(MouseEvent e) {

    }

    @Override
    public void mouseEntered(MouseEvent e) {

    }

    @Override
    public void mouseExited(MouseEvent e) {

    }

    @Override
    public void mouseDragged(MouseEvent e) {
        
        
        plotter_.cleanList();
        plotter_.plot("hugoniot", e.getX(),e.getY());
        paint(getGraphics());
        plotter_.update();
        
        
     
    }

    @Override
    public void mouseMoved(MouseEvent e) {
        
//          
//        plotter_.cleanList();
//        plotter_.plot("curva em NativeCanvas!!!", e.getX(),e.getY());
//        paint(getGraphics());
//        plotter_.update();
        
        

    }

    @Override
    public void componentResized(ComponentEvent e) {

    }

    @Override
    public void componentMoved(ComponentEvent e) {

    }

    @Override
    public void componentShown(ComponentEvent e) {

    }

    @Override
    public void componentHidden(ComponentEvent e) {

    }

}
