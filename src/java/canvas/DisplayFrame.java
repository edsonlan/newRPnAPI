/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package canvas;

import control.Plotter;
import java.awt.HeadlessException;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import javax.swing.JFrame;

/**
 *
 * @author edsonlan
 */
public class DisplayFrame extends JFrame implements WindowListener{
    
    
    private NativeCanvas canvas_;

    public DisplayFrame(String title,int frameId) throws HeadlessException {
        super(title);
        addWindowListener(this);
        canvas_= new NativeCanvas(frameId,new Plotter());
        getContentPane().add(canvas_);
        setSize(400,400);

    }

    public NativeCanvas getCanvas() {
        return canvas_;
    }
    
    
    @Override
    public void windowOpened(WindowEvent e) {

    }

    @Override
    public void windowClosing(WindowEvent e) {
      

    }

    @Override
    public void windowClosed(WindowEvent e) {
        dispose();
        System.exit(0);
    }

    @Override
    public void windowIconified(WindowEvent e) {

    }

    @Override
    public void windowDeiconified(WindowEvent e) {

    }

    @Override
    public void windowActivated(WindowEvent e) {

    }

    @Override
    public void windowDeactivated(WindowEvent e) {

    }

    
    
}
