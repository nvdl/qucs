/***************************************************************************
                              mouseactions.cpp
                             ------------------
    begin                : Thu Aug 28 2003
    copyright            : (C) 2003 by Michael Margraf
    email                : michael.margraf@alumni.tu-berlin.de
 ***************************************************************************/

/* Copyright (C) 2014 Guilherme Brondani Torri <guitorri@gmail.com>        */

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "qucs.h"
#include "main.h"
#include "node.h"
#include "schematic.h"
#include "mouseactions.h"
#include "module.h"
#include "components/component.h"
#include "components/spicedialog.h"
#include "components/spicefile.h"
#include "components/optimizedialog.h"
#include "components/componentdialog.h"
#include "components/vacomponent.h"
#include "components/sp_customsim.h"
#include "diagrams/diagramdialog.h"
#include "diagrams/markerdialog.h"
#include "diagrams/tabdiagram.h"
#include "diagrams/timingdiagram.h"
#include "dialogs/labeldialog.h"
#include "extsimkernels/customsimdialog.h"

#include <QTextStream>
#include <Q3PtrList>
#include <QMouseEvent>
#include <QClipboard>
#include <QApplication>
#include <QMessageBox>
#include <QMenu>
#include <QEvent>
#include <QAction>
#include <QLineEdit>
#include <QDebug>

#include <limits.h>
#include <stdlib.h>


#define  DOC_X_POS(x)  (int(float(x)/Doc->Scale) + Doc->ViewX1)
#define  DOC_Y_POS(y)  (int(float(y)/Doc->Scale) + Doc->ViewY1)
#define  DOC_X_FPOS    (float(Event->pos().x())/Doc->Scale + float(Doc->ViewX1))
#define  DOC_Y_FPOS    (float(Event->pos().y())/Doc->Scale + float(Doc->ViewY1))

#define  SCR_X_POS(x)  int(float(x - Doc->ViewX1) * Doc->Scale)
#define  SCR_Y_POS(y)  int(float(y - Doc->ViewY1) * Doc->Scale)

QAction *formerAction;   // remember action before drag n'drop etc.


MouseActions::MouseActions(QucsApp* App_)
{
  App = App_; // pointer to main app
  selElem  = 0;  // no component/diagram is selected
  isMoveEqual = false;  // mouse cursor move x and y the same way
  focusElement = 0; //element being interacted with mouse

  // ...............................................................
  // initialize menu appearing by right mouse button click on component
  ComponentMenu = new QMenu(QucsMain);
  focusMEvent   = new QMouseEvent(QEvent::MouseButtonPress, QPoint(0,0),
				  Qt::NoButton, Qt::NoButton);
}


MouseActions::~MouseActions()
{
  delete ComponentMenu;
  delete focusMEvent;
}

// -----------------------------------------------------------
void MouseActions::setPainter(Schematic *Doc)
{
  // contents to viewport transformation

  Doc->PostPaintEvent(_Translate,-Doc->contentsX(), -Doc->contentsY());
  Doc->PostPaintEvent(_Scale,Doc->Scale, Doc->Scale);
  Doc->PostPaintEvent(_Translate,-Doc->ViewX1, -Doc->ViewY1);
  Doc->PostPaintEvent(_DotLine);
  Doc->PostPaintEvent(_NotRop);

}

// -----------------------------------------------------------
bool MouseActions::pasteElements(Schematic *Doc)
{
  QClipboard *cb = QApplication::clipboard();   // get system clipboard
  QString s = cb->text(QClipboard::Clipboard);
  QTextStream stream(&s, QIODevice::ReadOnly);
  movingElements.clear();
  if(!Doc->paste(&stream, &movingElements)) return false;

  Element *pe;
  int xmax, xmin, ymax, ymin;
  xmin = ymin = INT_MAX;
  xmax = ymax = INT_MIN;
  // First, get the max and min coordinates of all selected elements.
  for(pe = movingElements.first(); pe != 0; pe = movingElements.next()) {
    if(pe->Type == isWire) {
      if(pe->x1 < xmin) xmin = pe->x1;
      if(pe->x2 > xmax) xmax = pe->x2;
      if(pe->y1 < ymin) ymin = pe->y1;
      if(pe->y2 > ymax) ymax = pe->y2;
    }
    else {
      if(pe->cx < xmin) xmin = pe->cx;
      if(pe->cx > xmax) xmax = pe->cx;
      if(pe->cy < ymin) ymin = pe->cy;
      if(pe->cy > ymax) ymax = pe->cy;
    }
  }

  xmin = -((xmax+xmin) >> 1);   // calculate midpoint
  ymin = -((ymax+ymin) >> 1);
  Doc->setOnGrid(xmin, ymin);

  // moving with mouse cursor in the midpoint
  for(pe = movingElements.first(); pe != 0; pe = movingElements.next())
    if(pe->Type & isLabel) {
      pe->cx += xmin;  pe->x1 += xmin;
      pe->cy += ymin;  pe->y1 += ymin;
    }
    else
      pe->setCenter(xmin, ymin, true);

  return true;
}

// -----------------------------------------------------------
void MouseActions::editLabel(Schematic *Doc, WireLabel *pl)
{
  LabelDialog *Dia = new LabelDialog(pl, Doc);
  int Result = Dia->exec();
  if(Result == 0) return;

  QString Name  = Dia->NodeName->text();
  QString Value = Dia->InitValue->text();
  delete Dia;

  if(Name.isEmpty() && Value.isEmpty()) { // if nothing entered, delete label
    pl->pOwner->Label = 0;   // delete name of wire
    delete pl;
  }
  else {
/*    Name.replace(' ', '_');	// label must not contain spaces
    while(Name.at(0) == '_') Name.remove(0,1);  // must not start with '_'
    if(Name.isEmpty()) return;
    if(Name == pl->Name) return;*/
    if(Result == 1) return;  // nothing changed

    int old_x2 = pl->x2;
    pl->setName(Name);   // set new name
    pl->initValue = Value;
    if(pl->cx > (pl->x1+(pl->x2>>1)))
      pl->x1 -= pl->x2 - old_x2; // don't change position due to text width
  }

  Doc->sizeOfAll(Doc->UsedX1, Doc->UsedY1, Doc->UsedX2, Doc->UsedY2);
  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
// Reinserts all elements (moved by the user) back into the schematic.
void MouseActions::endElementMoving(Schematic *Doc, Q3PtrList<Element> *movElements)
{
  Element *pe;
  for(pe = movElements->first(); pe!=0; pe = movElements->next()) {
//    pe->isSelected = false;  // deselect first (maybe afterwards pe == NULL)
    switch(pe->Type) {
      case isWire:
        if(pe->x1 == pe->x2)
          if(pe->y1 == pe->y2) {
            // Delete wires with zero length, but preserve label.
            if(((Wire*)pe)->Label) {
              Doc->insertNodeLabel((WireLabel*)((Wire*)pe)->Label);
              ((Wire*)pe)->Label = 0;
            }
            delete (Wire*)pe;
            break;
          }

	Doc->insertWire((Wire*)pe);
	break;
      case isDiagram:
	Doc->Diagrams->append((Diagram*)pe);
	break;
      case isPainting:
	Doc->Paintings->append((Painting*)pe);
	break;
      case isComponent:
      case isAnalogComponent:
      case isDigitalComponent:
	Doc->insertRawComponent((Component*)pe, false);
	break;
      case isMovingLabel:
      case isHMovingLabel:
      case isVMovingLabel:
	Doc->insertNodeLabel((WireLabel*)pe);
	break;
      case isMarker:
	((Marker*)pe)->pGraph->Markers.append((Marker*)pe);
	break;
    }
  }

  movElements->clear();
  if((MAx3 != 0) || (MAy3 != 0))  // moved or put at the same place ?
    Doc->setChanged(true, true);

  // enlarge viewarea if components lie outside the view
  Doc->sizeOfAll(Doc->UsedX1, Doc->UsedY1, Doc->UsedX2, Doc->UsedY2);
  Doc->enlargeView(Doc->UsedX1, Doc->UsedY1, Doc->UsedX2, Doc->UsedY2);
  Doc->viewport()->update();
  //drawn = false;
}

// -----------------------------------------------------------
// Moves elements in "movElements" by x/y
void MouseActions::moveElements(Q3PtrList<Element> *movElements, int x, int y)
{

  qDebug() << "moveElements";

  Wire *pw;
  Element *pe;

  for (pe = movElements->first(); pe != 0; pe = movElements->next()) {
    if (pe->Type == isWire) {
      pw = (Wire*) pe;   // connected wires are not moved completely

      if (((unsigned long)pw->Port1) > 3) {
        pw->x1 += x;  pw->y1 += y;

        if (pw->Label) {
          pw->Label->cx += x;  pw->Label->cy += y;
        }

      } else {
        if (long(pw->Port1) & 1) {
          pw->x1 += x;
        }

        if (long(pw->Port1) & 2) {
          pw->y1 += y;
        }
      }

      if (((unsigned long) pw->Port2) > 3) {
        pw->x2 += x;  pw->y2 += y;
      } else {
        if (long(pw->Port2) & 1) pw->x2 += x;
        if (long(pw->Port2) & 2) pw->y2 += y;
      }

      if (pw->Label) {      // root of node label must lie on wire
        if (pw->Label->cx < pw->x1) pw->Label->cx = pw->x1;
        if (pw->Label->cy < pw->y1) pw->Label->cy = pw->y1;
        if (pw->Label->cx > pw->x2) pw->Label->cx = pw->x2;
        if (pw->Label->cy > pw->y2) pw->Label->cy = pw->y2;
      }

    } else {
      pe->setCenter(x, y, true);
    }
  }
}
// ***********************************************************************
// **********                                                   **********
// **********       Functions for serving mouse moving          **********
// **********                                                   **********
// ***********************************************************************
void MouseActions::MMoveElement(Schematic *Doc, QMouseEvent *Event)
{
  qDebug() << "MMoveElement";

  if (selElem == 0) return;

  Doc->grabKeyboard(); // nvdl: For space key rotation

//  qDebug() << "MMoveElement got selElem";

  int x  = Event->pos().x();
  int y  = Event->pos().y();
  int fx = DOC_X_POS(x);
  int fy = DOC_Y_POS(y);
  int gx = fx;
  int gy = fy;
  Doc->setOnGrid(gx, gy);

  //QPainter painter(Doc->viewport());
  setPainter(Doc);

  if (selElem->Type == isPainting) {
    Doc->PostPaintEvent (_NotRop, 0,0,0,0);
    x -= Doc->contentsX();
    y -= Doc->contentsY();
    ((Painting*)selElem)->MouseMoving(Doc, x, y, gx, gy,
                                       Doc, x, y, drawn);
    drawn = true;
    Doc->viewport()->update();
    return;
  }  // of "if(isPainting)"


  // ********** it is a component or diagram
  /*if (drawn) selElem->paintScheme(Doc); // erase old scheme
  drawn = true;
  Doc->viewport()->repaint();*/

//  Component *comp = (Component*)selElem;
  //qDebug() << "desc" << comp->Description << "gx" << gx << "gy" << gy;

  selElem->setCenter(gx, gy);
  selElem->paintScheme(Doc); // paint scheme at new position

  //Doc->viewport()->repaint();
  Doc->viewport()->update();
}

// -----------------------------------------------------------
/**
 * @brief MouseActions::MMoveWire2 Paint wire as it is being drawn with mouse.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveWire2(Schematic *Doc, QMouseEvent *Event)
{

  //qDebug() << "MMoveWire2";

  if (Event != 0) { // If not just called for the wire orientation toggle update
    MAx2 = DOC_X_POS(Event->pos().x());
    MAy2 = DOC_Y_POS(Event->pos().y());
    Doc->setOnGrid(MAx2, MAy2);
  }

  drawWire(Doc, true); // Draw wire outline

  MCloseToNode(Doc, Event); // Draw cross if close to a node

  //QucsMain->MouseDoubleClickAction = &MouseActions::MDoubleClickWire2;
  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveWire1 Paint hair cross for "insert wire" mode
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveWire1(Schematic *Doc, QMouseEvent *Event)
{

  //qDebug() << "MMoveWire1";

  MAx2 = DOC_X_POS(Event->pos().x());
  MAy2 = DOC_Y_POS(Event->pos().y());

  Doc->setOnGrid(MAx2, MAy2);

  /*MAx2  = DOC_X_POS(Doc->viewport()->width());
  MAy2  = DOC_Y_POS(Doc->viewport()->height());
  Doc->PostPaintEvent (_Line, Doc->ViewX1, MAy3, MAx2, MAy3);
  Doc->PostPaintEvent (_Line, MAx3, Doc->ViewY1, MAx3, MAy2);*/
  
  if (!MCloseToNode(Doc, Event)) {
	  //nvdl: todo: Make the cross size and color user configurable
	  Doc->PostPaintEvent(_Line, MAx2 - 25, MAy2, MAx2 + 25, MAy2);
	  Doc->PostPaintEvent(_Line, MAx2, MAy2 - 25, MAx2, MAy2 + 25);
  }
  
  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveSelect Paints a selection rectangle or moves a selection.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveSelect(Schematic *Doc, QMouseEvent *Event)
{
  int x1, y1, width, height;
  //int xMin = 10000, xMax = 0, yMin = 10000, yMax = 0;
  //Element *pe;

  if (focusElement) {
    // print define value in hex, see element.h
    qDebug() << "MMoveSelect: focusElement->Type" <<  QString("0x%1").arg(focusElement->Type, 0, 16);

    /*for (pe=movingElements.first(); pe!=0; pe=movingElements.next()) {
        if (pe->x1 < xMin) xMin = pe->x1;
        if (pe->y1 < yMin) yMin = pe->y1;
        if (pe->x2 > xMax) xMax = pe->x2;
        if (pe->y2 > yMax) yMax = pe->y2;
    }

    width = xMax - xMin;
    height = yMax - yMin;

    Doc->PostPaintEvent(_Rect, x1, y1, width, height);*/

    MAx2 = DOC_X_POS(Event->pos().x()) - MAx1;
    MAy2 = DOC_Y_POS(Event->pos().y()) - MAy1;

    Doc->setOnGrid(MAx2, MAy2);

    if (isMoveEqual) { // Check if x and y movements should be equal
      if (abs(MAx2) > abs(MAy2)) {
        if (MAx2 < 0) MAx2 = -abs(MAy2); else MAx2 = abs(MAy2);
      } else {
        if (MAy2 < 0) MAy2 = -abs(MAx2); else MAy2 = abs(MAx2);
      }
    }

    Doc->PostPaintEvent (_Rect, MAx1, MAy1, MAx2, MAy2);

  } else {
    qDebug() << "MMoveSelect: Nothing under the mouse";

    MAx2 = DOC_X_POS(Event->pos().x());
    MAy2 = DOC_Y_POS(Event->pos().y());

    if (MAx1 <  MAx2) x1 = MAx1;
    else x1 = MAx2;

    if (MAy1 <  MAy2) y1 = MAy1;
    else y1 = MAy2;

    width = abs(MAx1 - MAx2);
    height = abs(MAy1 - MAy2);

    Doc->PostPaintEvent(_Rect, x1, y1, width, height);
  }

  //Doc->viewport()->repaint();

  //QucsMain->MouseMoveAction = &MouseActions::MMoveMoving;

  //qDebug() << "MMoveSelect " << "select area";

  /*MAx2 = DOC_X_POS(Event->pos().x()) - MAx1;
  MAy2 = DOC_Y_POS(Event->pos().y()) - MAy1;

  if(isMoveEqual) {    // x and y size must be equal ?
    if(abs(MAx2) > abs(MAy2)) {
      if(MAx2<0) MAx2 = -abs(MAy2); else MAx2 = abs(MAy2);
    }
    else { if(MAy2<0) MAy2 = -abs(MAx2); else MAy2 = abs(MAx2); }
  }

  Doc->PostPaintEvent (_Rect, MAx1, MAy1, MAx2, MAy2);*/
}

// -----------------------------------------------------------
void MouseActions::MMoveResizePainting(Schematic *Doc, QMouseEvent *Event)
{
  setPainter(Doc);

  MAx1 = DOC_X_POS(Event->pos().x());
  MAy1 = DOC_Y_POS(Event->pos().y());
  Doc->setOnGrid(MAx1, MAy1);
  ((Painting*)focusElement)->MouseResizeMoving(MAx1, MAy1, Doc);
}
/**
 * @brief MouseActions::MMoveMoving Moves components while the mouse button is kept pressed
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveMoving(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MMoveMoving";

  int centerX, centerY;

  setPainter(Doc);

  Doc->grabKeyboard();

  MAx2 = DOC_X_POS(Event->pos().x());
  MAy2 = DOC_Y_POS(Event->pos().y());

  Doc->setOnGrid(MAx2, MAy2);

  MAx3 = MAx1 = MAx2 - MAx1;
  MAy3 = MAy1 = MAy2 - MAy1;

  movingElements.clear();
  Doc->copySelectedElements(&movingElements);

  Wire *pw;
  // Changes the position of all moving elements by dx/dy

  for (Element *pe=movingElements.first(); pe!=0; pe=movingElements.next()) {

    if (pe->Type == isWire) {
      pw = (Wire*)pe;   // connecting wires are not moved completely

      if (((unsigned long)pw->Port1) > 3) {
        pw->x1 += MAx1;  pw->y1 += MAy1;
      } else {
        if (long(pw->Port1) & 1) {
          pw->x1 += MAx1;
        }

        if (long(pw->Port1) & 2) {
          pw->y1 += MAy1;
        }
      }

      if (((unsigned long)pw->Port2) > 3) {
        pw->x2 += MAx1;  pw->y2 += MAy1;
      } else {
        if (long(pw->Port2) & 1) pw->x2 += MAx1;
        if (long(pw->Port2) & 2) pw->y2 += MAy1;
      }

      if (pw->Label) {      // root of node label must lie on wire
        if (pw->Label->cx < pw->x1) pw->Label->cx = pw->x1;
        if (pw->Label->cy < pw->y1) pw->Label->cy = pw->y1;
        if (pw->Label->cx > pw->x2) pw->Label->cx = pw->x2;
        if (pw->Label->cy > pw->y2) pw->Label->cy = pw->y2;
      }

    } else {
      // nvdl: Align a component to the grid if off-grid
      pe->getCenter(centerX, centerY);
      Doc->setOnGrid(centerX, centerY);
      pe->setCenter(centerX, centerY, false);
      //pe->setCenter(MAx1, MAy1, true);
    }

    pe->paintScheme(Doc);
    qDebug() << "MMoveMoving: paintScheme";
  }

  drawn = true;
  MAx1 = MAx2;
  MAy1 = MAy2;
  QucsMain->MouseMoveAction = &MouseActions::MMoveMoving2;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseMoving;

  Doc->viewport()->repaint();
}

// -----------------------------------------------------------
// Moves components by keeping the mouse button pressed.
void MouseActions::MMoveMoving2(Schematic *Doc, QMouseEvent *Event)
{
  qDebug() << "MMoveMoving2";

  setPainter(Doc);

  MAx2 = DOC_X_POS(Event->pos().x());
  MAy2 = DOC_Y_POS(Event->pos().y());

  Doc->grabKeyboard();

  //if(drawn) // erase old scheme
    //for(pe = movingElements.first(); pe != 0; pe = movingElements.next())
      //pe->paintScheme(Doc);
//      if(pe->Type == isWire)  if(((Wire*)pe)->Label)
//        if(!((Wire*)pe)->Label->isSelected)
//          ((Wire*)pe)->Label->paintScheme(&painter);

  //drawn = true;

  // Use grid when CTRL key is not pressed
  if ((Event->state() & Qt::ControlModifier) == 0) {
    Doc->setOnGrid(MAx2, MAy2);
  }

  MAx1 = MAx2 - MAx1;
  MAy1 = MAy2 - MAy1;
  MAx3 += MAx1;  MAy3 += MAy1; // keep track of the complete movement

  moveElements(&movingElements, MAx1, MAy1);  // moves elements by MAx1/MAy1
  paintElementsScheme(Doc);

  //if(pe->Type == isWire)  if(((Wire*)pe)->Label)
  //if(!((Wire*)pe)->Label->isSelected)
  //  ((Wire*)pe)->Label->paintScheme(&painter);

  MAx1 = MAx2;
  MAy1 = MAy2;

  Doc->viewport()->repaint();
  //Doc->viewport()->update();
}

// -----------------------------------------------------------
/**
 * @brief MouseActions::MMovePaste Moves components after paste from clipboard.
 * @param Doc
 * @param Event
 */
void MouseActions::MMovePaste(Schematic *Doc, QMouseEvent *Event)
{
  MAx1 = DOC_X_POS(Event->pos().x());
  MAy1 = DOC_Y_POS(Event->pos().y());
  moveElements(Doc,MAx1,MAy1);
  paintElementsScheme(Doc);

  drawn = true;
  QucsMain->MouseMoveAction = &MouseActions::MMoveMoving2;
  QucsMain->MouseReleaseAction = &MouseActions::MReleasePaste;

  Doc->viewport()->update();
}

// -----------------------------------------------------------
// Moves scroll bar of diagram (e.g. tabular) according the mouse cursor.
void MouseActions::MMoveScrollBar(Schematic *Doc, QMouseEvent *Event)
{
  TabDiagram *d = (TabDiagram*)focusElement;
  int x = DOC_X_POS(Event->pos().x());
  int y = DOC_Y_POS(Event->pos().y());

  if(d->scrollTo(MAx2, x - MAx1, y - MAy1)) {
    Doc->setChanged(true, true, 'm'); // 'm' = only the first time

// FIXME #warning QPainter p(Doc->viewport());
    // FIXME #warning ViewPainter Painter;
    // FIXME #warning Painter.init(&p, Doc->Scale, -Doc->ViewX1, -Doc->ViewY1,
// FIXME #warning                  Doc->contentsX(), Doc->contentsY());
// FIXME #warning     Painter.fillRect(d->cx-d->x1, d->cy-d->y2, d->x2+d->x1, d->y2+d->y1,
// FIXME #warning                      QucsSettings.BGColor);
// FIXME #warning     d->paint(&Painter);
  }
}

// -----------------------------------------------------------
/**
* @brief MouseActions::MMoveDelete
*   Paints a cross under the mouse cursor to show the delete mode.
* @param Doc Schematic document
* @param Event
*/
void MouseActions::MMoveDelete(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  // cannot draw on the viewport, it is displaced by the size of dock and toolbar
  Doc->PostPaintEvent (_Line, MAx3-15, MAy3-15, MAx3+15, MAy3+15,0,0,false);
  Doc->PostPaintEvent (_Line, MAx3-15, MAy3+15, MAx3+15, MAy3-15,0,0,false);

  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveLabel Paints a label above the mouse cursor for "set wire label".
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveLabel(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  // paint marker
  Doc->PostPaintEvent (_Line, MAx3, MAy3, MAx3+10, MAy3-10);
  Doc->PostPaintEvent (_Line, MAx3+10, MAy3-10, MAx3+20, MAy3-10);
  Doc->PostPaintEvent (_Line, MAx3+10, MAy3-10, MAx3+10, MAy3-17);

  // paint A
  Doc->PostPaintEvent (_Line, MAx3+12, MAy3-12, MAx3+15, MAy3-23);
  Doc->PostPaintEvent (_Line, MAx3+14, MAy3-17, MAx3+17, MAy3-17);
  Doc->PostPaintEvent (_Line, MAx3+19, MAy3-12, MAx3+16, MAy3-23);

  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveMarker Paints a triangle above the mouse for "set marker on graph"
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveMarker(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3, MAy3-2, MAx3-8, MAy3-10);
  Doc->PostPaintEvent (_Line, MAx3+1, MAy3-3, MAx3+8, MAy3-10);
  Doc->PostPaintEvent (_Line, MAx3-7, MAy3-10, MAx3+7, MAy3-10);

  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveMirrorX Paints rounded "mirror about y axis" mouse cursor
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveMirrorY(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3-11, MAy3-4, MAx3-9, MAy3-9);
  Doc->PostPaintEvent (_Line, MAx3-11, MAy3-3, MAx3-6, MAy3-3);
  Doc->PostPaintEvent (_Line, MAx3+11, MAy3-4, MAx3+9, MAy3-9);
  Doc->PostPaintEvent (_Line, MAx3+11, MAy3-3, MAx3+6, MAy3-3);
  Doc->PostPaintEvent (_Arc, MAx3-10, MAy3-8, 21, 10, 16*20, 16*140,false);

  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveMirrorX Paints rounded "mirror about x axis" mouse cursor
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveMirrorX(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3-4, MAy3-11, MAx3-9, MAy3-9);
  Doc->PostPaintEvent (_Line, MAx3-3, MAy3-11, MAx3-3, MAy3-6);
  Doc->PostPaintEvent (_Line, MAx3-4, MAy3+11, MAx3-9, MAy3+9);
  Doc->PostPaintEvent (_Line, MAx3-3, MAy3+11, MAx3-3, MAy3+6);
  Doc->PostPaintEvent (_Arc, MAx3-8, MAy3-10, 10, 21, 16*110, 16*140,false);

  Doc->viewport()->update();
}

/**
 * @brief MouseActions::MMoveMirrorX Paints "rotate" mouse cursor
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveRotate(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3-6, MAy3+8, MAx3-6, MAy3+1);
  Doc->PostPaintEvent (_Line, MAx3-7, MAy3+8, MAx3-12, MAy3+8);
  Doc->PostPaintEvent (_Arc, MAx3-10, MAy3-10, 21, 21, -16*20, 16*240,false);

  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveActivate Paints a crossed box mouse cursor to "(de)activate" components.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveActivate(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MMoveActivate";

  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Rect, MAx3, MAy3-9, 14, 10);
  Doc->PostPaintEvent (_Line, MAx3, MAy3-9, MAx3+13, MAy3);
  Doc->PostPaintEvent (_Line, MAx3, MAy3, MAx3+13, MAy3-9);

  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveOnGrid Paints a grid besides the mouse cursor, put "on grid" mode.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveOnGrid(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3+10, MAy3+ 3, MAx3+25, MAy3+3);
  Doc->PostPaintEvent (_Line, MAx3+10, MAy3+ 7, MAx3+25, MAy3+7);
  Doc->PostPaintEvent (_Line, MAx3+10, MAy3+11, MAx3+25, MAy3+11);
  Doc->PostPaintEvent (_Line, MAx3+13, MAy3, MAx3+13, MAy3+15);
  Doc->PostPaintEvent (_Line, MAx3+17, MAy3, MAx3+17, MAy3+15);
  Doc->PostPaintEvent (_Line, MAx3+21, MAy3, MAx3+21, MAy3+15);

  Doc->viewport()->update();
}

/**
 * @brief MouseActions::MMoveFreely Called when there is no active operation.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveFreely(Schematic *Doc, QMouseEvent *Event) {

	//qDebug() << "MMoveFreely";

	// nvdl: todo: Temporary fix; find out when drawing operations complete.
  // Escape key press?
  //Doc->releaseKeyboard();

  MAx2 = DOC_X_POS(Event->pos().x());
  MAy2 = DOC_Y_POS(Event->pos().y());

  Doc->setOnGrid(MAx2, MAy2);

  if (QucsSettings.NodeWiring) {
    if (MCloseToNode(Doc, Event)) {
      //qDebug() << "Close to node";
    }
  }

	Doc->viewport()->update();
}

/**
 * @brief MouseActions::MCloseToNode Checks for a near-by node and draws a cross if a node is found.
 * @param Doc
 * @param Event
 */
bool MouseActions::MCloseToNode(Schematic *Doc, QMouseEvent *Event) {

  //qDebug() << "MCloseToNode";

	Node *pn;
	bool nodeFound = false;
	int snapDistance = QucsSettings.grid1Spacing / 2;

	if (Doc->Nodes == NULL) {
	  return false;
	} else if (Doc->Nodes->isEmpty()) {
    return false;
  }

	for (pn = Doc->Nodes->first(); pn != 0; pn = Doc->Nodes->next()) {
		  if (abs(pn->cx - MAx2) <= snapDistance && abs(pn->cy - MAy2) <= snapDistance) {
		    MAx2 = pn->cx; // Allow off-grid node snap
		    MAy2 = pn->cy;
	    //if (pn->cx == MAx2 && pn->cy == MAy2) {
			  nodeFound = true;
			  break;
		  }
	}

	if (nodeFound) {
		//qDebug() << "MCloseToNode: Found a node";
		Doc->PostPaintEvent (_Line, MAx2 - 20, MAy2 - 20, MAx2 + 20, MAy2 + 20);
		Doc->PostPaintEvent (_Line, MAx2 - 20, MAy2 + 20, MAx2 + 20, MAy2 - 20);
		return true;
	} else {
		return false;
	}
}

/**
 * @brief MouseActions::MMoveMoveTextB Paints mouse symbol for "move component text" mode.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveMoveTextB(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3+14, MAy3   , MAx3+16, MAy3);
  Doc->PostPaintEvent (_Line, MAx3+23, MAy3   , MAx3+25, MAy3);
  Doc->PostPaintEvent (_Line, MAx3+13, MAy3   , MAx3+13, MAy3+ 3);
  Doc->PostPaintEvent (_Line, MAx3+13, MAy3+ 7, MAx3+13, MAy3+10);
  Doc->PostPaintEvent (_Line, MAx3+14, MAy3+10, MAx3+16, MAy3+10);
  Doc->PostPaintEvent (_Line, MAx3+23, MAy3+10, MAx3+25, MAy3+10);
  Doc->PostPaintEvent (_Line, MAx3+26, MAy3   , MAx3+26, MAy3+ 3);
  Doc->PostPaintEvent (_Line, MAx3+26, MAy3+ 7, MAx3+26, MAy3+10);

  //Doc->viewport()->update();
}


/**
 * @brief MouseActions::MMoveMoveText Paint rectangle around component text being mouse moved
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveMoveText(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MMoveMoveText";

  int newX = DOC_X_POS(Event->pos().x());
  int newY = DOC_Y_POS(Event->pos().y());

  MAx1 += newX - MAx3;
  MAy1 += newY - MAy3;
  MAx3 = newX;
  MAy3 = newY;

  //Doc->PostPaintEvent (_Rect, newX - 10, newY - 10, 20, 20);

  ((Component*)focusElement)->tx = newX - ((Component*)focusElement)->cx - 10;
  ((Component*)focusElement)->ty = newY - ((Component*)focusElement)->cy - 10;

  Doc->viewport()->update();

  //drawn = false;
  Doc->setChanged(true, true);
}


/**
 * @brief MouseActions::MMoveZoomIn Paints symbol beside the mouse to show the "Zoom in" modus.
 * @param Doc
 * @param Event
 */
void MouseActions::MMoveZoomIn(Schematic *Doc, QMouseEvent *Event)
{
  MAx3  = DOC_X_POS(Event->pos().x());
  MAy3  = DOC_Y_POS(Event->pos().y());

  Doc->PostPaintEvent (_Line, MAx3+14, MAy3   , MAx3+22, MAy3);
  Doc->PostPaintEvent (_Line, MAx3+18, MAy3-4 , MAx3+18, MAy3+4);
  Doc->PostPaintEvent (_Ellipse, MAx3+12, MAy3-6, 13, 13,0,0,false);

  //Doc->viewport()->update();
}

/**
 * @brief MouseActions::MMovePanning Pans the schematic view.
 * @param Doc
 * @param Event
 */
void MouseActions::MMovePanning(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MMovePanning";

  float xShift;
  float yShift;

  MAx2 = Event->pos().x() - Doc->contentsX();
  MAy2 = Event->pos().y() - Doc->contentsY();

  xShift = (MAx2 - MAx1);
  yShift = (MAy2 - MAy1);

  /*qDebug() << "MMovePanning: ViewX1: " << Doc->ViewX1;
  qDebug() << "MMovePanning: ViewX2: " << Doc->ViewX2;
  qDebug() << "MMovePanning: ViewY1: " << Doc->ViewY1;
  qDebug() << "MMovePanning: ViewY2: " << Doc->ViewY2;

  qDebug() << "MMovePanning: MAx1: " << MAx1;
  qDebug() << "MMovePanning: MAy1: " << MAy1;
  qDebug() << "MMovePanning: MAx2: " << MAx2;
  qDebug() << "MMovePanning: MAy2: " << MAy2;

  qDebug() << "MMovePanning: xShift: " << xShift;
  qDebug() << "MMovePanning: yShift: " << yShift;

  qDebug() << "MMovePanning: Doc->contentsX(): " << Doc->contentsX();
  qDebug() << "MMovePanning: Doc->contentsX(): " << Doc->contentsY();*/

  MAx1 = MAx2;
  MAy1 = MAy2;

  Doc->scrollBy(-xShift, 0);
  Doc->scrollBy(0, -yShift);

  panningDone = true;

  //Doc->viewport()->update();
  Doc->viewport()->repaint();
}

/**
 * @brief MouseActions::MReleasePanning Stops panning of the schematic view.
 * @param Doc
 * @param Event
 */
void MouseActions::MReleasePanning(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleasePanning";

  QucsMain->MouseMoveAction = &MouseActions::MMoveFreely;
  QucsMain->MousePressAction = &MouseActions::MPressSelect;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect;
  QucsMain->MouseDoubleClickAction = &MouseActions::MDoubleClickSelect;

  panningDone = false;

  //Doc->releaseKeyboard();
}

// ************************************************************************
// **********                                                    **********
// **********    Functions for serving mouse button clicking     **********
// **********                                                    **********
// ************************************************************************

// Is called from several MousePress functions to show right button menu.
void MouseActions::rightPressMenu(Schematic *Doc, QMouseEvent *Event, float fX, float fY)
{

  int x1, y1;

  qDebug() << "rightPressMenu";

  if (panningDone && QucsMain->MouseMoveAction == &MouseActions::MMovePanning) {
    qDebug() << "rightPressMenu: Panning stopped";
    MReleasePanning(Doc, Event); // Switch mouse call-backs to normal
    return;

  } else if (QucsMain->MouseMoveAction == &MouseActions::MMoveMoving2) { // Rotation
    x1 = DOC_X_POS(Event->pos().x());
    y1 = DOC_Y_POS(Event->pos().y());
    rotateElements(Doc, x1, y1);
    paintElementsScheme(Doc);
    //Doc->viewport()->update();

    return;

  } else if (QucsMain->MousePressAction == &MouseActions::MPressElement) {
    return;

  } else if (QucsMain->MousePressAction == &MouseActions::MPressWire2) {
    return;

  } else if (QucsMain->MouseMoveAction == &MouseActions::MMovePanning) {
    MReleasePanning(Doc, Event); // No panning happened, switch mouse call-backs to normal
  }

  qDebug() << "rightPressMenu: Opening menu";

  MAx1 = int(fX);
  MAy1 = int(fY);
  focusElement = Doc->selectElement(fX, fY, false);

  if(focusElement)  // remove special function (4 least significant bits)
    focusElement->Type &= isSpecialMask;


  // define menu
  ComponentMenu->clear();
  while(true) {
    if(focusElement) {
      focusElement->isSelected = true;
      ComponentMenu->insertItem(
         QObject::tr("Edit Properties"), QucsMain, SLOT(slotEditElement()));

      if((focusElement->Type & isComponent) == 0) break;
    }
    else {
/// \todo "exchange like this"
      //ComponentMenu->addAction(QucsMain->symEdit);
      //to QucsMain->symEdit->addTo(ComponentMenu);
      // see http://qt-project.org/doc/qt-4.8/qaction-qt3.html#addTo
      QucsMain->symEdit->addTo(ComponentMenu);
      QucsMain->fileSettings->addTo(ComponentMenu);
    }
    if(!QucsMain->moveText->isOn())
      QucsMain->moveText->addTo(ComponentMenu);
    break;
  }

  // nvdl: todo: Add timeouts or other break checks
  while(true) {
    if(focusElement)
      if(focusElement->Type == isGraph) break;
    if(!QucsMain->onGrid->isOn())
      QucsMain->onGrid->addTo(ComponentMenu);
    QucsMain->editCopy->addTo(ComponentMenu);
    if(!QucsMain->editPaste->isOn())
      QucsMain->editPaste->addTo(ComponentMenu);
    break;
  }

  while (true) {
    if (focusElement) {
      if (focusElement->Type == isDiagram) {
        ComponentMenu->insertItem(QObject::tr("Export as image"), QucsMain,
            SLOT(slotSaveDiagramToGraphicsFile()));
      }
    }
    break;
  }

  if(!QucsMain->editDelete->isOn())
    QucsMain->editDelete->addTo(ComponentMenu);
  if(focusElement) if(focusElement->Type == isMarker) {
    ComponentMenu->insertSeparator();
    QString s = QObject::tr("power matching");
    if( ((Marker*)focusElement)->pGraph->Var == "Sopt" )
      s = QObject::tr("noise matching");
    ComponentMenu->insertItem(s, QucsMain, SLOT(slotPowerMatching()));
    if( ((Marker*)focusElement)->pGraph->Var.left(2) == "S[" )
      ComponentMenu->insertItem(QObject::tr("2-port matching"), QucsMain,
                                SLOT(slot2PortMatching()));
  }
  do {
    if(focusElement) {
      if(focusElement->Type == isDiagram) break;
      if(focusElement->Type == isGraph) {
        QucsMain->graph2csv->addTo(ComponentMenu);
        break;
      }
    }
    ComponentMenu->insertSeparator();
    if(focusElement) if(focusElement->Type & isComponent)
      if(!QucsMain->editActivate->isOn())
        QucsMain->editActivate->addTo(ComponentMenu);
    if(!QucsMain->editRotate->isOn())
      QucsMain->editRotate->addTo(ComponentMenu);
    if(!QucsMain->editMirror->isOn())
      QucsMain->editMirror->addTo(ComponentMenu);
    if(!QucsMain->editMirrorY->isOn())
      QucsMain->editMirrorY->addTo(ComponentMenu);

    // right-click menu to go into hierarchy
    if(focusElement) {
      if(focusElement->Type & isComponent)
	if(((Component*)focusElement)->Model == "Sub")
	  if(!QucsMain->intoH->isOn())
	    QucsMain->intoH->addTo(ComponentMenu);
    }
    // right-click menu to pop out of hierarchy
    if(!focusElement)
      if(!QucsMain->popH->isOn())
	QucsMain->popH->addTo(ComponentMenu);
  } while(false);

  *focusMEvent = *Event;  // remember event for "edit component" action
  ComponentMenu->popup(Event->globalPos());
  //drawn = false;

  Doc->viewport()->update();
}

// -----------------------------------------------------------
void MouseActions::MPressLabel(Schematic *Doc, QMouseEvent*, float fX, float fY)
{

  qDebug() << "MPressLabel";

  int x = int(fX), y = int(fY);
  Wire *pw = 0;
  WireLabel *pl=0;
  Node *pn = Doc->selectedNode(x, y);
  if(!pn) {
    pw = Doc->selectedWire(x, y);
    if(!pw) return;
  }

  QString Name, Value;
  Element *pe=0;
  // is wire line already labeled ?
  if(pw) pe = Doc->getWireLabel(pw->Port1);
  else pe = Doc->getWireLabel(pn);
  if(pe) {
    if(pe->Type & isComponent) {
      QMessageBox::information(0, QObject::tr("Info"),
                 QObject::tr("The ground potential cannot be labeled!"));
      return;
    }
    pl = ((Conductor*)pe)->Label;
  }

  LabelDialog *Dia = new LabelDialog(pl, Doc);
  if(Dia->exec() == 0) return;

  Name  = Dia->NodeName->text();
  Value = Dia->InitValue->text();
  delete Dia;

  if(Name.isEmpty() && Value.isEmpty() ) { // if nothing entered, delete name
    if(pe) {
      if(((Conductor*)pe)->Label)
        delete ((Conductor*)pe)->Label; // delete old name
      ((Conductor*)pe)->Label = 0;
    }
    else {
      if(pw) pw->setName("", "");   // delete name of wire
      else pn->setName("", "");
    }
  }
  else {
/*    Name.replace(' ', '_');	// label must not contain spaces
    while(Name.at(0) == '_') Name.remove(0,1);  // must not start with '_'
    if(Name.isEmpty()) return;
*/
    if(pe) {
      if(((Conductor*)pe)->Label)
        delete ((Conductor*)pe)->Label; // delete old name
      ((Conductor*)pe)->Label = 0;
    }

    int xl = x+30;
    int yl = y-30;
    Doc->setOnGrid(xl, yl);
    // set new name
    if(pw) pw->setName(Name, Value, x-pw->x1 + y-pw->y1, xl, yl);
    else pn->setName(Name, Value, xl, yl);
  }

  Doc->sizeOfAll(Doc->UsedX1, Doc->UsedY1, Doc->UsedX2, Doc->UsedY2);
  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
void MouseActions::MPressSelect(Schematic *Doc, QMouseEvent *Event, float fX, float fY)
{
  MAx1 = int(fX);
  MAy1 = int(fY);

  if (Event->button() == Qt::MiddleButton) {

    qDebug() << "MPressSelect: Middle mouse press";
    return;

  } else if (Event->button() == Qt::RightButton) {

    qDebug() << "MPressSelect: Right mouse press";

    if (QucsMain->MouseMoveAction == 0 || QucsMain->MouseMoveAction == &MouseActions::MMoveFreely) {
      QucsMain->MouseMoveAction = &MouseActions::MMovePanning;
      QucsMain->MouseReleaseAction = &MouseActions::MReleasePanning;
      QucsMain->MousePressAction = 0;
      QucsMain->MouseDoubleClickAction = 0;

      //Doc->grabKeyboard();

      // For panning
      MAx1 = Event->pos().x() - Doc->contentsX();
      MAy1 = Event->pos().y() - Doc->contentsY();
    }

    panningDone = false;

    return;
  }

  QucsMain->MouseReleaseAction = 0;
  QucsMain->MouseMoveAction = 0;
  QucsMain->MousePressAction = 0;
  QucsMain->MouseDoubleClickAction = 0;

  bool Ctrl;
  if (Event->state() & Qt::ControlModifier) Ctrl = true;
  else Ctrl = false;

  int No = 0;

  focusElement = Doc->selectElement(fX, fY, Ctrl, &No);
  isMoveEqual = false;   // moving not neccessarily square

  if (focusElement)
    // print define value in hex, see element.h
    qDebug() << "MPressSelect: focusElement->Type" <<  QString("0x%1").arg(focusElement->Type, 0, 16);
  else
    qDebug() << "MPressSelect: Nothing under the mouse";

  drawn = true;

  if (focusElement) {

    switch(focusElement->Type)
    {
      case isPainting:
        QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect;
        QucsMain->MouseMoveAction = &MouseActions::MMoveMoving;
        //QucsMain->MousePressAction = 0;
        //QucsMain->MouseDoubleClickAction = 0;
        Doc->grabKeyboard();  // no keyboard inputs during move actions
        // Update matching wire label highlighting
        Doc->highlightWireLabels ();
        break;

      case isPaintingResize:  // resize painting ?
        focusElement->Type = isPainting;
        QucsMain->MouseReleaseAction = &MouseActions::MReleaseResizePainting;
        QucsMain->MouseMoveAction = &MouseActions::MMoveResizePainting;
        //QucsMain->MousePressAction = 0;
        //QucsMain->MouseDoubleClickAction = 0;
        Doc->grabKeyboard();  // no keyboard inputs during move actions
        // Update matching wire label highlighting
        Doc->highlightWireLabels ();
        break;

      case isDiagramResize:  // resize diagram ?
        if(((Diagram*)focusElement)->Name.left(4) != "Rect")
          if(((Diagram*)focusElement)->Name.at(0) != 'T')
            if(((Diagram*)focusElement)->Name != "Curve")
              isMoveEqual = true;  // diagram must be square

        focusElement->Type = isDiagram;
        MAx1 = focusElement->cx;
        MAx2 = focusElement->x2;

        if (((Diagram*)focusElement)->State & 1) {
          MAx1 += MAx2;
          MAx2 *= -1;
        }

        MAy1 =  focusElement->cy;
        MAy2 = -focusElement->y2;

        if (((Diagram*)focusElement)->State & 2) {
          MAy1 += MAy2;
          MAy2 *= -1;
        }

        QucsMain->MouseReleaseAction = &MouseActions::MReleaseResizeDiagram;
        QucsMain->MouseMoveAction = &MouseActions::MMoveSelect;
        //QucsMain->MousePressAction = 0;
        //QucsMain->MouseDoubleClickAction = 0;

        Doc->grabKeyboard(); // no keyboard inputs during move actions
        // Update matching wire label highlighting
        Doc->highlightWireLabels ();

        break;

      case isDiagramHScroll:  // scroll in tabular ?
        MAy1 = MAx1;

      case isDiagramVScroll:
        focusElement->Type = isDiagram;

        No = ((TabDiagram*)focusElement)->scroll(MAy1);

        switch(No)
        {
          case 1:
            Doc->setChanged(true, true, 'm'); // 'm' = only the first time
            break;
          case 2:  // move scroll bar with mouse cursor
            QucsMain->MouseMoveAction = &MouseActions::MMoveScrollBar;
            //QucsMain->MousePressAction = 0;
            //QucsMain->MouseDoubleClickAction = 0;
            Doc->grabKeyboard();  // no keyboard inputs during move actions

            // Remember inital scroll bar position.
            MAx2 = int(((TabDiagram*)focusElement)->xAxis.limit_min);
            // Update matching wire label highlighting
            Doc->highlightWireLabels ();
            break;
        }
        // Update matching wire label highlighting
        Doc->highlightWireLabels();
        Doc->viewport()->update();
        //drawn = false;
        break;

      case isComponentText:  // property text of component ?

        if (Ctrl) {
          focusElement->Type &= (~isComponentText) | isComponent;
          MAx3 = No;
          QucsMain->slotApplyCompText();
          // Update matching wire label highlighting
          Doc->highlightWireLabels();
        } else {
          MPressMoveText(Doc, Event, fX, fY);
        }

        break;

      case isNode:
        if (QucsSettings.NodeWiring) {
          formerAction = QucsMain->select; // to restore action afterwards
          QucsMain->activeAction = QucsMain->insWire;

          QucsMain->select->blockSignals(true);
          QucsMain->select->setChecked(false);
          QucsMain->select->blockSignals(false);

          QucsMain->insWire->blockSignals(true);
          QucsMain->insWire->setChecked(true);
          QucsMain->insWire->blockSignals(false);

          // Update matching wire label highlighting
          Doc->highlightWireLabels ();

          toggleWireOrientation = false;
          MPressWire1(Doc, Event, fX, fY);
        }

        break;

      case isWire:
        QucsMain->MouseMoveAction = &MouseActions::MMoveMoving; //MMoveWire2;
        QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect;
        //QucsMain->MousePressAction = &MouseActions::MPressWire2;
        break;

      default:
        qDebug() << "MPressSelect: Unknown element";
        QucsMain->MouseMoveAction = &MouseActions::MMoveMoving;
        QucsMain->MouseReleaseAction = &MouseActions::MReleaseMoving;
    }
  }

  //if (!drawn) {
    //QucsMain->MousePressAction = 0;
    //QucsMain->MouseDoubleClickAction = 0;
    //Doc->grabKeyboard();  // no keyboard inputs during move actions
    //Doc->viewport()->update();
  //}

  if (focusElement == 0) {
    //MAx2 = 0;  // if not clicking on an element => open a rectangle
    //MAy2 = 0;

    if (QucsMain->MouseReleaseAction == 0) {
      QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect2;
    }

    // If there is no element, start a selection rectangle
    if (QucsMain->MouseMoveAction == 0) {
      QucsMain->MouseMoveAction = &MouseActions::MMoveSelect;
    }

    ////drawn = false;

  } else {
    // Element can be moved
    if (!Ctrl) { // If Ctrl key is not pressed, only select the clicked one
      if (!focusElement->isSelected) // Don't move selected elements if clicked
        Doc->deselectElements(focusElement); // Deselect all elements except "focusElement"

      focusElement->isSelected = true;
    }

    Doc->setOnGrid(MAx1, MAy1);
    //QucsMain->MouseMoveAction = &MouseActions::MMoveMoving;
  }

  // Update matching wire label highlighting
  Doc->highlightWireLabels();

  Doc->viewport()->update();
}

// -----------------------------------------------------------
void MouseActions::MPressDelete(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  Element *pe = Doc->selectElement(fX, fY, false);
  if(pe)
  {
    pe->isSelected = true;
    Doc->deleteElements();

    Doc->sizeOfAll(Doc->UsedX1, Doc->UsedY1, Doc->UsedX2, Doc->UsedY2);
    Doc->viewport()->update();
    //drawn = false;
  }
}

// -----------------------------------------------------------
void MouseActions::MPressActivate(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  MAx1 = int(fX);
  MAy1 = int(fY);
  if(!Doc->activateSpecifiedComponent(MAx1, MAy1)) {
//    if(Event->button() != Qt::LeftButton) return;
    MAx2 = 0;  // if not clicking on a component => open a rectangle
    MAy2 = 0;
    QucsMain->MousePressAction = 0;
    QucsMain->MouseReleaseAction = &MouseActions::MReleaseActivate;
    QucsMain->MouseMoveAction = &MouseActions::MMoveSelect;
  }

  //drawn = false;
  Doc->viewport()->update();
}

// -----------------------------------------------------------
void MouseActions::MPressMirrorX(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  // no use in mirroring wires or diagrams
  Component *c = Doc->selectedComponent(int(fX), int(fY));
  if(c) {
    if(c->Ports.count() < 1) return;  // only mirror components with ports
    c->mirrorX();
    Doc->setCompPorts(c);
  }
  else {
    Painting *p = Doc->selectedPainting(fX, fY);
    if(p == 0) return;
    p->mirrorX();
  }

  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
void MouseActions::MPressMirrorY(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  // no use in mirroring wires or diagrams
  Component *c = Doc->selectedComponent(int(fX), int(fY));
  if(c) {
    if(c->Ports.count() < 1) return;  // only mirror components with ports
    c->mirrorY();
    Doc->setCompPorts(c);
  }
  else {
    Painting *p = Doc->selectedPainting(fX, fY);
    if(p == 0) return;
    p->mirrorY();
  }

  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
void MouseActions::MPressRotate(Schematic *Doc, QMouseEvent*, float fX, float fY)
{

  qDebug() << "MPressRotate";

  int centerX, centerY;
  Element *e = Doc->selectElement(int(fX), int(fY), false);
  if (e == 0) return;
  e->Type &= isSpecialMask;  // remove special functions

  WireLabel *pl;
  int x1, y1, x2, y2;
//  e->isSelected = false;
  switch(e->Type) {
    case isComponent:
    case isAnalogComponent:
    case isDigitalComponent:
      if(((Component*)e)->Ports.count() < 1)
        break;  // do not rotate components without ports
      ((Component*)e)->rotate();
      Doc->setCompPorts((Component*)e);
      // enlarge viewarea if component lies outside the view
      ((Component*)e)->entireBounds(x1,y1,x2,y2, Doc->textCorr());
      Doc->enlargeView(x1, y1, x2, y2);
      break;

    case isWire:
      pl = ((Wire*)e)->Label;
      ((Wire*)e)->Label = 0;    // prevent label to be deleted
      Doc->Wires->setAutoDelete(false);
      Doc->deleteWire((Wire*)e);
      ((Wire*)e)->Label = pl;
      ((Wire*)e)->rotate();
      Doc->setOnGrid(e->x1, e->y1);
      Doc->setOnGrid(e->x2, e->y2);
      if(pl)  Doc->setOnGrid(pl->cx, pl->cy);
      Doc->insertWire((Wire*)e);
      Doc->Wires->setAutoDelete(true);
      if (Doc->Wires->containsRef ((Wire*)e))
        Doc->enlargeView(e->x1, e->y1, e->x2, e->y2);
      break;

    case isPainting:
      ((Painting*)e)->rotate();
      // enlarge viewarea if component lies outside the view
      ((Painting*)e)->Bounding(x1,y1,x2,y2);
      Doc->enlargeView(x1, y1, x2, y2);
      break;

    default:
      qDebug() << "MPressRotate: Unknown element";
      return;
  }

  // nvdl: todo: Align center to grid (needs fixing as nodes get detached from the component)
  /*switch (e->Type) {
  case isComponent:
  case isAnalogComponent:
  case isDigitalComponent:
    e->getCenter(centerX, centerY);
    Doc->setOnGrid(centerX, centerY);
    e->setCenter(centerX, centerY, false);
    break;

  default:
    break;
  }*/

  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
// insert component, diagram, painting into schematic ?!
void MouseActions::MPressElement(Schematic *Doc, QMouseEvent *Event, float, float)
{
  qDebug() << "MPressElement";

  if (selElem == 0) return;
  //QPainter painter(Doc->viewport());
  //setPainter(Doc, &painter);

  int x1, y1, x2, y2, rot;

  if (selElem->Type & isComponent) {
    Component *Comp = (Component*)selElem;
    //qDebug() << "+-+ got to switch:" << Comp->Name;
    QString entryName = Comp->Name;

    switch (Event->button()) {

    case Qt::LeftButton:
      // left mouse button inserts component into the schematic
      // give the component a pointer to the schematic it's a
      // part of
      Comp->setSchematic(Doc);
      Comp->textSize(x1, y1);
      Doc->insertComponent(Comp);
      Comp->textSize(x2, y2);
      if(Comp->tx < Comp->x1) Comp->tx -= x2 - x1;

      // Note: insertCopmponents does increment  name1 -> name2
      //qDebug() << "  +-+ got to insert:" << Comp->Name;

      // enlarge viewarea if component lies outside the view
      Comp->entireBounds(x1,y1,x2,y2, Doc->textCorr());
      Doc->enlargeView(x1, y1, x2, y2);

      //drawn = false;
      //Doc->viewport()->update();
      Doc->setChanged(true, true);
      rot = Comp->rotated;

      // handle static and dynamic components
      //QucsApp::CompChoose;

      if (Module::vaComponents.contains(entryName)){
        QString filename = Module::vaComponents[entryName];
        //qDebug() << "   ===+ recast";
        Comp = dynamic_cast<vacomponent*>(Comp)->newOne(filename); //va component
        qDebug() << "   => recast = Comp;" << Comp->Name << "filename: " << filename;
      } else {
        Comp = Comp->newOne(); // static component is used, so create a new one
      }

      rot -= Comp->rotated;
      rot &= 3;

      // keep last rotation for single component
      while (rot--) Comp->rotate();

      break;

    case Qt::RightButton:  // Right mouse button rotates the component

      if (Comp->Ports.count() == 0)
        break;  // do not rotate components without ports

      //Comp->paintScheme(Doc); // erase old component scheme
      Comp->rotate();
      //Doc->viewport()->repaint();
      Comp->paintScheme(Doc); // paint new component scheme

      break;

    default: ;   // avoids compiler warnings
    }

    // qDebug() << "   => selElem = Comp;" << Comp->Name;
    // comp it geting empty
    selElem = Comp;
    return;

  } else if (selElem->Type == isDiagram) {
    if (Event->button() != Qt::LeftButton) return;

    Diagram *Diag = (Diagram*) selElem;
    QFileInfo Info(Doc->DocName);

    // dialog is Qt::WDestructiveClose !!!
    DiagramDialog *dia =
    new DiagramDialog(Diag, Doc);

    if (dia->exec() == QDialog::Rejected) {  // don't insert if dialog canceled
      Doc->viewport()->update();
      //drawn = false;
      return;
    }

    Doc->Diagrams->append(Diag);
    Doc->enlargeView(Diag->cx, Diag->cy-Diag->y2, Diag->cx+Diag->x2, Diag->cy);
    Doc->setChanged(true, true);   // document has been changed

    //Doc->viewport()->repaint();

    Diag = Diag->newOne(); // the component is used, so create a new one
    Diag->paintScheme(Doc);
    selElem = Diag;
    Doc->viewport()->update();
    return;
  }  // of "if(isDiagram)"

  // ***********  it is a painting !!!
  if (((Painting*)selElem)->MousePressing()) {
    Doc->Paintings->append((Painting*)selElem);
    ((Painting*)selElem)->Bounding(x1,y1,x2,y2);
    //Doc->enlargeView(x1, y1, x2, y2);
    selElem = ((Painting*)selElem)->newOne();

    //Doc->viewport()->update();
    Doc->setChanged(true, true);

    MMoveElement(Doc, Event);  // needed before next mouse pressing
    //drawn = false;
  }

  Doc->viewport()->update();
}

/**
 * @brief MouseActions::MReleaseElement Is called when there is an element to be placed and mouse receives release event
 * @param Doc
 * @param Event
  */
void MouseActions::MReleaseElement(Schematic *Doc, QMouseEvent *Event) {

  Doc->releaseKeyboard();
  Doc->viewport()->update();
}

/**
 * @brief MouseActions::MPressWire1 Is called if starting point of wire is pressed
 * @param Doc
 * @param Event
 * @param fX
 * @param fY
 */
void MouseActions::MPressWire1(Schematic *Doc, QMouseEvent *Event, float fX, float fY)
{

  qDebug() << "MPressWire1";

  toggleWireOrientation = false;

  Doc->grabKeyboard(); // To capture keyboard shortcuts

  MAx3 = int(fX);
  MAy3 = int(fY);
  Doc->setOnGrid(MAx3, MAy3);

  // No previous wire to follow the orientation
  prevWireX1 = 0;
  prevWireY1 = 0;
  prevWireX2 = 0;
  prevWireY2 = 0;

  formerAction = 0; // keep wire action active after first wire finished
  QucsMain->MouseMoveAction = &MouseActions::MMoveWire2;
  QucsMain->MousePressAction = &MouseActions::MPressWire2;
  QucsMain->MouseReleaseAction = 0;
  QucsMain->MouseDoubleClickAction = &MouseActions::MDoubleClickWire2;

  // Double-click action is set in "MMoveWire2" to not initiate it
  // during "Wire1" actions.
  Doc->viewport()->update();
}


/**
 * @brief MouseActions::MPressWire2 Is called if ending point of wire is pressed
 * @param Doc
 * @param Event
 * @param fX
 * @param fY
 */
void MouseActions::MPressWire2(Schematic *Doc, QMouseEvent *Event, float fX, float fY)
{

  qDebug() << "MPressWire2";

  MAx2 = int(fX);
  MAy2 = int(fY);
  Doc->setOnGrid(MAx2, MAy2);

  //Doc->releaseKeyboard();

  switch (Event->button()) {

  case Qt::LeftButton:
    drawWire(Doc, false); // Draw real wire(s)
    break;

   /// \todo document right mouse button changes the wire corner
  case Qt::RightButton:

    //drawWire(Doc, true); // Draw wire outline
    //Doc->viewport()->repaint();

    /*MAx2 = int(fX);
    MAy2 = int(fY);
    Doc->setOnGrid(MAx2, MAy2);*/

    //MAx1 ^= 1; // Change the direction of wire corner
    toggleWireOrientation = not toggleWireOrientation;
    drawWire(Doc, true); // Draw the new outline

    break;

  default: ;    // avoids compiler warnings
  }

  Doc->viewport()->update();
}

// -----------------------------------------------------------
/**
 * @brief MouseActions::drawWire Draws an actual wire or its outline
 * @param Doc
 * @param Event
 * @param fX
 * @param fY
 */
void MouseActions::drawWire(Schematic* Doc, bool outline) {

  int dx, dy, xIntersect, yIntersect, xLen, yLen;
  int set1 = 0, set2 = 0;
  bool drawVerHor, wireFinished;
  bool prevWireExists, prevWireVertical;

  dx = abs(prevWireX2 - prevWireX1);
  dy = abs(prevWireY2 - prevWireY1);

  if (dx == 0 && dy == 0) { // No previous wire, follow the mouse pointer
    dx = abs(MAx3 - MAx2);
    dy = abs(MAy3 - MAy2);
    //qDebug() << "drawWire: No previous wire to follow";
    prevWireExists = false;
  } else {
    prevWireExists = true;
    if (dy > dx) {
      prevWireVertical = true;
    } else {
      prevWireVertical = false;
    }
  }

  if (dy > dx) {
    drawVerHor = true;
    //qDebug() << "drawWire: Vertical wire";
  } else {
    drawVerHor = false;
    //qDebug() << "drawWire: Horizontal wire";
  }

  if (toggleWireOrientation) drawVerHor = not drawVerHor; // Toggle order

  // Force order to avoid "go-back" on the same wire
  if (prevWireExists) {
    if (prevWireVertical) {
      if (prevWireY2 < prevWireY1 && MAy2 > prevWireY2) drawVerHor = false;
      if (prevWireY2 > prevWireY1 && MAy2 < prevWireY2) drawVerHor = false;
    } else {
      if (prevWireX2 < prevWireX1 && MAx2 > prevWireX2) drawVerHor = true;
      if (prevWireX2 > prevWireX1 && MAx2 < prevWireX2) drawVerHor = true;
    }
  }

  if (drawVerHor) {
    xIntersect = MAx3;
    yIntersect = MAy2;
  } else {
    xIntersect = MAx2;
    yIntersect = MAy3;
  }

  if (outline) { // Draw wire(s) outlines(s)
    Doc->PostPaintEvent (_Line, MAx3, MAy3, xIntersect, yIntersect);
    Doc->PostPaintEvent (_Line, xIntersect, yIntersect, MAx2, MAy2);

  } else { // Draw actual wire(s)
    xLen = abs(MAx3 - MAx2);
    yLen = abs(MAy3 - MAy2);

    if (drawVerHor) {
      if (xLen > 0) {
        set1 = Doc->insertWire(new Wire(xIntersect, yIntersect, MAx2, MAy2));
      }

      if (yLen > 0) {
        set2 = set1;
        set1 = Doc->insertWire(new Wire(MAx3, MAy3, xIntersect, yIntersect));
      }

    } else {
      if (xLen > 0) {
        set1 = Doc->insertWire(new Wire(MAx3, MAy3, xIntersect, yIntersect));
      }

      if (yLen > 0) {
        set2 = set1;
        set1 = Doc->insertWire(new Wire(xIntersect, yIntersect, MAx2, MAy2));
      }
    }

    prevWireX1 = xIntersect;
    prevWireY1 = yIntersect;
    prevWireX2 = MAx2;
    prevWireY2 = MAy2;

    /*qDebug() << "prevWireX1:" << prevWireX1;
    qDebug() << "prevWireY1:" << prevWireY1;
    qDebug() << "prevWireX2:" << prevWireX2;
    qDebug() << "prevWireY2:" << prevWireY2;*/

    QString s1 = QString("0x%1").arg(set1, 0, 16);
    QString s2 = QString("0x%1").arg(set2, 0, 16);

    qDebug() << "set1: " << s1;
    qDebug() << "set2: " << s2;

    wireFinished = false;

    // nvdl: todo: More testing as "set1" and "set2" values are not well understood
    if (set2 == 0) { // Only a vertical or horizontal wire was drawn
      if ((set1 & 0x3) == 0x3) {
        wireFinished = true;
      }
    } else {
      if ((set1 & 0x3) == 0x3 && (set2 & 0x2) == 0x2) {
        wireFinished = true;
      }
    }

    //if ((set1 & 0x3) == 0x3 || (set2 & 0x3) == 0x3) { // if(set1 & 2) {
    if (!wireFinished) {
      // if last port is connected, then...
      if (formerAction) { // nvdl: todo: Make more sense of it
        // Restore old action
        QucsMain->select->setChecked(true);

        qDebug() << "drawWire: formerAction";

      } else {
        // Start a new wire
        QucsMain->MouseMoveAction = &MouseActions::MMoveWire2;
        QucsMain->MousePressAction = &MouseActions::MPressWire2;
        QucsMain->MouseDoubleClickAction = 0;

        qDebug() << "drawWire: Continuing same wire";
      }

    } else {
      qDebug() << "drawWire: Wire finished";

      QucsMain->MouseMoveAction = &MouseActions::MMoveWire1;
      QucsMain->MousePressAction = &MouseActions::MPressWire1;
      QucsMain->MouseDoubleClickAction = 0;
    }

    ////drawn = false; // nvdl: todo: Legacy?
    if (set1 | set2) Doc->setChanged(true, true);
    MAx3 = MAx2;
    MAy3 = MAy2;

    toggleWireOrientation = false;
  }
}
// -----------------------------------------------------------
// Is called for setting a marker on a diagram's graph
void MouseActions::MPressMarker(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  MAx1 = int(fX);
  MAy1 = int(fY);
  Marker *pm = Doc->setMarker(MAx1, MAy1);

  if(pm) {
    int x0 = pm->Diag->cx;
    int y0 = pm->Diag->cy;
    Doc->enlargeView(x0+pm->x1, y0-pm->y1-pm->y2, x0+pm->x1+pm->x2, y0-pm->y1);
  }
  Doc->viewport()->update();
  //drawn = false;
}

// -----------------------------------------------------------
void MouseActions::MPressOnGrid(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  Element *pe = Doc->selectElement(fX, fY, false);
  if(pe)
  {
    pe->Type &= isSpecialMask;  // remove special functions (4 lowest bits)

    // onGrid is toggle action -> no other element can be selected
    pe->isSelected = true;
    Doc->elementsOnGrid();

    Doc->sizeOfAll(Doc->UsedX1, Doc->UsedY1, Doc->UsedX2, Doc->UsedY2);
    // Update matching wire label highlighting
    Doc->highlightWireLabels ();
    Doc->viewport()->update();
    //drawn = false;
  }

}

// -----------------------------------------------------------
void MouseActions::MPressMoveText(Schematic *Doc, QMouseEvent*, float fX, float fY)
{

  qDebug() << "MPressMoveText";

  MAx1 = int(fX);
  MAy1 = int(fY);
  focusElement = Doc->selectCompText(MAx1, MAy1, MAx2, MAy2);

  if (focusElement) {
    MAx3 = MAx1;
    MAy3 = MAy1;
    MAx1 = ((Component*)focusElement)->cx + ((Component*)focusElement)->tx;
    MAy1 = ((Component*)focusElement)->cy + ((Component*)focusElement)->ty;
    //Doc->viewport()->update();
    //drawn = false;
    QucsMain->MouseMoveAction = &MouseActions::MMoveMoveText;
    QucsMain->MouseReleaseAction = &MouseActions::MReleaseMoveText;
    //Doc->grabKeyboard();  // no keyboard inputs during move actions
  }
}

// -----------------------------------------------------------
void MouseActions::MPressZoomIn(Schematic *Doc, QMouseEvent*, float fX, float fY)
{
  qDebug() << "MPressZoomIn: Zoom into box";

  MAx1 = int(fX);
  MAy1 = int(fY);
  MAx2 = 0;  // rectangle size
  MAy2 = 0;

  QucsMain->MouseMoveAction = &MouseActions::MMoveSelect;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseZoomIn;
  Doc->grabKeyboard();  // no keyboard inputs during move actions
  Doc->viewport()->update();
  //drawn = false;
}


// ***********************************************************************
// **********                                                   **********
// **********    Functions for serving mouse button releasing   **********
// **********                                                   **********
// ***********************************************************************
void MouseActions::MReleaseSelect(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseSelect";

  bool ctrl;
  if(Event->state() & Qt::ControlModifier) ctrl = true;
  else ctrl = false;

  if(!ctrl) Doc->deselectElements(focusElement);

  if (focusElement && Event->button() == Qt::LeftButton) {
    if (focusElement->Type == isWire) {
      Doc->selectWireLine(focusElement, ((Wire*)focusElement)->Port1, ctrl);
      Doc->selectWireLine(focusElement, ((Wire*)focusElement)->Port2, ctrl);
    }
  }

  defaultState();
  Doc->releaseKeyboard();  // allow keyboard inputs again

  Doc->highlightWireLabels ();
  //drawn = false;

  Doc->viewport()->update();
}

// -----------------------------------------------------------
// Is called after the rectangle for selection is released.
void MouseActions::MReleaseSelect2(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseSelect2";

  Doc->releaseKeyboard(); // allow keyboard inputs again

  if(Event->button() != Qt::LeftButton) return;

  bool Ctrl;
  if (Event->state() & Qt::ControlModifier) Ctrl = true;
  else Ctrl = false;

  MAx2 = DOC_X_POS(Event->pos().x());
  MAy2 = DOC_Y_POS(Event->pos().y());

  // Select all elements within the rectangle
  Doc->selectElements(MAx1, MAy1, MAx2, MAy2, Ctrl);

  defaultState();

  Doc->highlightWireLabels ();
  //drawn = false;

  Doc->viewport()->update();
}

// -----------------------------------------------------------
void MouseActions::MReleaseActivate(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseActivate";

  if(Event->button() != Qt::LeftButton) return;

  // activates all components within the rectangle
  Doc->activateCompsWithinRect(MAx1, MAy1, MAx1+MAx2, MAy1+MAy2);

  QucsMain->MouseMoveAction = &MouseActions::MMoveActivate;
  QucsMain->MousePressAction = &MouseActions::MPressActivate;
  QucsMain->MouseReleaseAction = 0;
  QucsMain->MouseDoubleClickAction = 0;

  Doc->highlightWireLabels ();
  Doc->viewport()->update();
  //drawn = false;
}

// -----------------------------------------------------------
// Is called after moving selected elements.
void MouseActions::MReleaseMoving(Schematic *Doc, QMouseEvent*)
{

  qDebug() << "MReleaseMoving";

  // Allow all mouse buttons, because for others than the left one,
  // a menu has already created.
  endElementMoving(Doc, &movingElements);

  Doc->releaseKeyboard();  // allow keyboard inputs again
  defaultState();
}

// -----------------------------------------------------------
void MouseActions::MReleaseResizeDiagram(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseResizeDiagram";

  if(Event->button() != Qt::LeftButton) return;

  MAx3  = focusElement->cx;
  MAy3  = focusElement->cy;
  if(MAx2 < 0) {    // resize diagram
    if(MAx2 > -10) MAx2 = -10;   // not smaller than 10 pixels
    focusElement->x2 = -MAx2;
    focusElement->cx = MAx1+MAx2;
  }
  else {
    if(MAx2 < 10) MAx2 = 10;
    focusElement->x2 = MAx2;
    focusElement->cx = MAx1;
  }
  if(MAy2 < 0) {
    if(MAy2 > -10) MAy2 = -10;
    focusElement->y2 = -MAy2;
    focusElement->cy = MAy1;
  }
  else {
    if(MAy2 < 10) MAy2 = 10;
    focusElement->y2 = MAy2;
    focusElement->cy = MAy1+MAy2;
  }
  MAx3 -= focusElement->cx;
  MAy3 -= focusElement->cy;

  Diagram *pd = (Diagram*)focusElement;
  pd->updateGraphData();
  foreach(Graph *pg, pd->Graphs)
    foreach(Marker *pm, pg->Markers) {
      pm->x1 += MAx3;      // correct changes due to move of diagram corner
      pm->y1 += MAy3;
    }

  int x1, x2, y1, y2;
  pd->Bounding(x1, x2, y1, y2);
  Doc->enlargeView(x1, x2, y1, y2);

  QucsMain->MouseMoveAction = &MouseActions::MMoveFreely;
  QucsMain->MousePressAction = &MouseActions::MPressSelect;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect;
  QucsMain->MouseDoubleClickAction = &MouseActions::MDoubleClickSelect;
  //Doc->releaseKeyboard();  // allow keyboard inputs again

  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
void MouseActions::MReleaseResizePainting(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseResizePainting";

  if(Event->button() != Qt::LeftButton) return;

  QucsMain->MouseMoveAction = &MouseActions::MMoveFreely;
  QucsMain->MousePressAction = &MouseActions::MPressSelect;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect;
  QucsMain->MouseDoubleClickAction = &MouseActions::MDoubleClickSelect;
  //Doc->releaseKeyboard();  // allow keyboard inputs again

  Doc->viewport()->update();
  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
/**
 * @brief MouseActions::paintElementsScheme Draws the outline of all selected elements
 * @param Doc
 */
void MouseActions::paintElementsScheme(Schematic *Doc) {

  Element *pe;
  int x1, y1, width, height;
  int xMin = INT_MAX, xMax = INT_MIN, yMin = INT_MAX, yMax = INT_MIN;

  for (pe = movingElements.first(); pe != 0; pe = movingElements.next()) {
    if (pe->Type == isWire) {
      // nvdl: todo: Temporary fix until "paintScheme" of wire is fixed
      Doc->PostPaintEvent(_Line, pe->x1, pe->y1, pe->x2, pe->y2, 0, 0, false);
    } else {
      pe->paintScheme(Doc);
    }

    // nvdl: todo: Skip elements that are causing issues
    if (pe->Type & 0x4000 || pe->Type & 0x8000) {
      continue;
    }

    if (pe->cx + pe->x1 < xMin) xMin = pe->cx + pe->x1;
    if (pe->cy + pe->y1 < yMin) yMin = pe->cy + pe->y1;
    if (pe->cx + pe->x2 > xMax) xMax = pe->cx + pe->x2;
    if (pe->cy + pe->y2 > yMax) yMax = pe->cy + pe->y2;
  }

  /*for (pe = movingElements.first(); pe != 0; pe = movingElements.next()) {
      qDebug() << "xMin:" << xMin;
      qDebug() << "yMin:" << yMin;
      qDebug() << "xMax:" << xMax;
      qDebug() << "yMax:" << yMax;

      qDebug() << "pe-cx:" << pe->cx;
      qDebug() << "pe-cy:" << pe->cy;
      qDebug() << "pe-x1:" << pe->x1;
      qDebug() << "pe-y1:" << pe->y1;
      qDebug() << "pe-x2:" << pe->x2;
      qDebug() << "pe-y2:" << pe->y2;
      qDebug() << "pe-Type:" << pe->Type << endl;
  }*/

  xMin -= 25; yMin -= 25; xMax += 25; yMax += 25;

  width = xMax - xMin;
  height = yMax - yMin;

  Doc->PostPaintEvent(_Rect, xMin, yMin, width, height);
}
// -----------------------------------------------------------
void MouseActions::moveElements(Schematic *Doc, int& x1, int& y1)
{
  Element *pe;
  Doc->setOnGrid(x1, y1);

  for(pe=movingElements.first(); pe!=0; pe=movingElements.next()) {
    if(pe->Type & isLabel) {
      pe->cx += x1;  pe->x1 += x1;
      pe->cy += y1;  pe->y1 += y1;
    }
    else
      pe->setCenter(x1, y1, true);
  }
}

// -----------------------------------------------------------
void MouseActions::rotateElements(Schematic *Doc, int& x1, int& y1)
{

  qDebug() << "rotateElements";

  int x2, y2;
  int centerX, centerY;
  Element *pe;
  Doc->setOnGrid(x1, y1);

  // nvdl: todo: Updated the values (MAx1, MAy1) so that the component does not go off the grid.
  // It can happen if there is select->rotate->move situation.
  // Still there are issues.
  MAx1 = x1;
  MAy1 = y1;

  for (pe = movingElements.first(); pe != 0; pe = movingElements.next()) {
    switch (pe->Type) {
    case isComponent:
    case isAnalogComponent:
    case isDigitalComponent:
      ((Component*)pe)->rotate(); // rotate !before! setting the center
      x2 = x1 - pe->cx;
      pe->setCenter(pe->cy - y1 + x1, x2 + y1);
      break;

    case isWire:
      x2     = pe->x1;
      pe->x1 = pe->y1 - y1 + x1;
      pe->y1 = x1 - x2 + y1;
      x2     = pe->x2;
      pe->x2 = pe->y2 - y1 + x1;
      pe->y2 = x1 - x2 + y1;
      break;

    case isPainting:
      ((Painting*)pe)->rotate(); // rotate !before! setting the center
      ((Painting*)pe)->getCenter(x2, y2);
      pe->setCenter(y2 - y1 + x1, x1 - x2 + y1);
      break;

    default:
      x2 = x1 - pe->cx;   // if diagram -> only rotate cx/cy
      pe->setCenter(pe->cy - y1 + x1, x2 + y1);
      break;
    }

    // nvdl: Align center to grid
    switch (pe->Type) {
    case isComponent:
    case isAnalogComponent:
    case isDigitalComponent:
      pe->getCenter(centerX, centerY);
      Doc->setOnGrid(centerX, centerY);
      pe->setCenter(centerX, centerY, false);
      break;

    default:
      break;
    }

    /*if (pe->Type != isWire) {
      pe->getCenter(centerX, centerY);
      Doc->setOnGrid(centerX, centerY);
      pe->setCenter(centerX, centerY, false);
    }*/
  }
}

// -----------------------------------------------------------
void MouseActions::MReleasePaste(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleasePaste";

  int x1, y1, x2, y2, rot;
  QFileInfo Info(Doc->DocName);
  //QPainter painter(Doc->viewport());

  Element *pe;
  switch(Event->button()) {
  case Qt::LeftButton :
    // insert all moved elements into document
    for(pe = movingElements.first(); pe!=0; pe = movingElements.next()) {
      pe->isSelected = false;
      switch(pe->Type) {
	case isWire:
	  if(pe->x1 == pe->x2) if(pe->y1 == pe->y2)  break;
	  Doc->insertWire((Wire*)pe);
	  if (Doc->Wires->containsRef ((Wire*)pe))
	    Doc->enlargeView(pe->x1, pe->y1, pe->x2, pe->y2);
	  else pe = NULL;
	  break;
	case isDiagram:
	  Doc->Diagrams->append((Diagram*)pe);
	  ((Diagram*)pe)->loadGraphData(Info.dirPath() + QDir::separator() +
					Doc->DataSet);
	  Doc->enlargeView(pe->cx, pe->cy-pe->y2, pe->cx+pe->x2, pe->cy);
	  break;
	case isPainting:
	  Doc->Paintings->append((Painting*)pe);
	  ((Painting*)pe)->Bounding(x1,y1,x2,y2);
	  Doc->enlargeView(x1, y1, x2, y2);
	  break;
	case isMovingLabel:
	  pe->Type = isNodeLabel;
	  Doc->placeNodeLabel((WireLabel*)pe);
	  break;
	case isComponent:
	case isAnalogComponent:
	case isDigitalComponent:
	  Doc->insertComponent((Component*)pe);
	  ((Component*)pe)->entireBounds(x1,y1,x2,y2, Doc->textCorr());
	  Doc->enlargeView(x1, y1, x2, y2);
	  break;
      }
    }

    pasteElements(Doc);
    // keep rotation sticky for pasted elements
    rot = movingRotated;
    x1 = y1 = 0;

    while(rot--)
      rotateElements(Doc,x1,y1);

    //Doc->viewport()->repaint();
    paintElementsScheme(Doc);

    QucsMain->MouseMoveAction = &MouseActions::MMovePaste;
    QucsMain->MousePressAction = 0;
    QucsMain->MouseReleaseAction = 0;
    QucsMain->MouseDoubleClickAction = 0;

    //drawn = false;
    Doc->viewport()->update();
    Doc->setChanged(true, true);
    break;

  // ............................................................
  case Qt::RightButton :  // right button rotates the elements
    //setPainter(Doc, &painter);

    /*if(drawn) // erase old scheme // nvdl: todo: Possibly not needed (legacy)
      paintElementsScheme(Doc);

    drawn = true;*/

    x1 = DOC_X_POS(Event->pos().x());
    y1 = DOC_Y_POS(Event->pos().y());

    //Doc->viewport()->repaint();
    rotateElements(Doc,x1,y1);
    paintElementsScheme(Doc);
    // save rotation
    movingRotated++;
    movingRotated &= 3;
    break;

  default: ;    // avoids compiler warnings
  }

  //Doc->viewport()->update();
}

// -----------------------------------------------------------
void MouseActions::MReleaseMoveText(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseMoveText";

  if (Event->button() != Qt::LeftButton) return;

  defaultState();

  //QucsMain->MouseMoveAction = &MouseActions::MMoveMoveTextB;
  //QucsMain->MouseReleaseAction = 0;
  //Doc->releaseKeyboard();  // allow keyboard inputs again

  //((Component*)focusElement)->tx = MAx1 - ((Component*)focusElement)->cx;
  //((Component*)focusElement)->ty = MAy1 - ((Component*)focusElement)->cy;

  Doc->viewport()->update();

  //drawn = false;
  Doc->setChanged(true, true);
}

// -----------------------------------------------------------
void MouseActions::MReleaseZoomIn(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MReleaseZoomIn";

  if(Event->button() != Qt::LeftButton) return;

  MAx1 = Event->pos().x();
  MAy1 = Event->pos().y();
  float DX = float(MAx2);
  float DY = float(MAy2);

  float initialScale = Doc->Scale;
  float scale = 1;
  float xShift = 0;
  float yShift = 0;

  if((Doc->Scale * DX) < 6.0) {
    // a simple click zooms by constant factor
    scale = Doc->zoom(1.5)/initialScale;

    xShift = scale * Event->pos().x();
    yShift = scale * Event->pos().y();

  } else {
    float xScale = float(Doc->visibleWidth())  / abs(DX);
    float yScale = float(Doc->visibleHeight()) / abs(DY);
    scale = qMin(xScale, yScale)/initialScale;
    scale = Doc->zoom(scale)/initialScale;

    xShift = scale * (MAx1 - 0.5*DX);
    yShift = scale * (MAy1 - 0.5*DY);
  }

  xShift -= (0.5*Doc->visibleWidth() + Doc->contentsX());
  yShift -= (0.5*Doc->visibleHeight() + Doc->contentsY());

  Doc->scrollBy(xShift, yShift);

  QucsMain->MouseMoveAction = &MouseActions::MMoveZoomIn;
  QucsMain->MouseReleaseAction = 0;
  //Doc->releaseKeyboard();  // allow keyboard inputs again
}


// ***********************************************************************
// **********                                                   **********
// **********    Functions for mouse button double clicking     **********
// **********                                                   **********
// ***********************************************************************
void MouseActions::editElement(Schematic *Doc, QMouseEvent *Event)
{
//    qDebug() << "+double click, editElement";

  qDebug() << "editElement";

  if(focusElement == 0) return;

//  qDebug() << "+focusElement->Type" << focusElement->Type;

  Graph *pg;
  Component *c;
  Diagram *dia;
  DiagramDialog *ddia;
  MarkerDialog *mdia;
  int x1, y1, x2, y2;

  QFileInfo Info(Doc->DocName);
  float fX = DOC_X_FPOS, fY = DOC_Y_FPOS;

  switch(focusElement->Type) {
    case isComponent:
    case isAnalogComponent:
    case isDigitalComponent:
         c = (Component*)focusElement;
//         qDebug() << "cast focusElement into" << c->Name;
         if(c->Model == "GND") return;

         if (c->Model == ".CUSTOMSIM") {
             CustomSimDialog *sd = new CustomSimDialog((SpiceCustomSim*)c, Doc);
             if(sd->exec() != 1) break;   // dialog is WDestructiveClose
         } else if(c->Model == "SPICE") {
           SpiceDialog *sd = new SpiceDialog(App, (SpiceFile*)c, Doc);
           if(sd->exec() != 1) break;   // dialog is WDestructiveClose
         }
         else if(c->Model == ".Opt") {
           OptimizeDialog *od = new OptimizeDialog((Optimize_Sim*)c, Doc);
           if(od->exec() != 1) break;   // dialog is WDestructiveClose
         }
         else {
           ComponentDialog * cd = new ComponentDialog(c, Doc);
           if(cd->exec() != 1) break;   // dialog is WDestructiveClose

           Doc->Components->findRef(c);
           Doc->Components->take();
           Doc->setComponentNumber(c); // for ports/power sources
           Doc->Components->append(c);
         }

         Doc->setChanged(true, true);
         c->entireBounds(x1,y1,x2,y2, Doc->textCorr());
         Doc->enlargeView(x1,y1,x2,y2);
         break;

    case isDiagram :
         dia = (Diagram*)focusElement;

         if (dia->Name.at(0) == 'T') { // don't open dialog on scrollbar
           if (dia->Name == "Time") {
             if (dia->cy < int(fY)) {
               if (((TimingDiagram*)focusElement)->scroll(MAx1)) {
                 Doc->setChanged(true, true, 'm'); // 'm' = only the first time
               }
               break;
             }
           } else {
             if (dia->cx > int(fX)) {
               if (((TabDiagram*)focusElement)->scroll(MAy1)) {
                 Doc->setChanged(true, true, 'm'); // 'm' = only the first time
               }
               break;
             }
           }
         }

         ddia = new DiagramDialog(dia, Doc);
         if (ddia->exec() != QDialog::Rejected) { // is WDestructiveClose
           Doc->setChanged(true, true);
         }

         dia->Bounding(x1, x2, y1, y2);
         Doc->enlargeView(x1, x2, y1, y2);
         break;

    case isGraph :
      pg = (Graph*) focusElement;
      // searching diagram for this graph
      for (dia = Doc->Diagrams->last(); dia != 0; dia = Doc->Diagrams->prev()) {
        if (dia->Graphs.indexOf(pg) >= 0) {
          break;
        }
      }

      if(!dia) break;

      ddia = new DiagramDialog(dia, Doc, pg);

      if (ddia->exec() != QDialog::Rejected) {  // is WDestructiveClose
        Doc->setChanged(true, true);
      }

      break;

    case isWire:
      defaultState();
      MPressLabel(Doc, Event, fX, fY);
      break;

    case isNodeLabel:
    case isHWireLabel:
    case isVWireLabel:
      editLabel(Doc, (WireLabel*)focusElement);
      break;

    case isPainting:
      if ( ((Painting*)focusElement)->Dialog() ) {
        Doc->setChanged(true, true);
      }

      break;

    case isMarker:
      mdia = new MarkerDialog((Marker*)focusElement, Doc);
      if (mdia->exec() > 1) {
        Doc->setChanged(true, true);
      }

      break;

    default:
      break;
  }

  // Very strange: Now an open VHDL editor gets all the keyboard input !?!
  // I don't know why it only happens here, nor am I sure whether it only
  // happens here. Anyway, I hope the best and give the focus back to the
  // current document.
  Doc->setFocus();

  Doc->viewport()->update();
  //drawn = false;
}

// -----------------------------------------------------------
void MouseActions::MDoubleClickSelect(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MDoubleClickSelect";

  //Doc->releaseKeyboard();  // allow keyboard inputs again
  QucsMain->editText->setHidden(true);
  editElement(Doc, Event);

  defaultState();
}


/**
 * @brief MouseActions::MDoubleClickWire2  Double click terminates wire insertion.
 * @param Doc
 * @param Event
 */
void MouseActions::MDoubleClickWire2(Schematic *Doc, QMouseEvent *Event)
{

  qDebug() << "MDoubleClickWire2";

  MPressWire2(Doc, Event, DOC_X_FPOS, DOC_Y_FPOS);

  if(formerAction)
    QucsMain->select->setChecked(true);  // restore old action
  else {
    QucsMain->MouseMoveAction = &MouseActions::MMoveWire1;
    QucsMain->MousePressAction = &MouseActions::MPressWire1;
    QucsMain->MouseDoubleClickAction = 0;
  }
}
//=================================================================================================
/**
 * @brief MouseActions::keyPressEvent Keyboard key press event sent by the schematic.
 * @param Doc
 * @param Event
 */
void MouseActions::keyPressEvent(Schematic *doc, QKeyEvent *event) {

  qDebug() << "keyPressEvent";
}
//=================================================================================================
/**
 * @brief MouseActions::keyReleaseEvent Keyboard key release event sent by the schematic.
 * @param Doc
 * @param Event
 */
void MouseActions::keyReleaseEvent(Schematic *doc, QKeyEvent *event) {

  qDebug() << "keyReleaseEvent: event->key()" << event->key();

  switch (event->key()) {

  // nvdl: todo: Keys will come from the shortcut manager
  case Qt::Key_Space:

    if (QucsMain->MouseMoveAction == &MouseActions::MMoveMoving2) { // Rotation
      //x1 = DOC_X_POS(Event->pos().x());
      //y1 = DOC_Y_POS(Event->pos().y());
      rotateElements(doc, MAx1, MAy1);
      paintElementsScheme(doc);

    } else if (QucsMain->MouseMoveAction == &MouseActions::MMoveWire2) { // Wire orientation toggle
      toggleWireOrientation = not toggleWireOrientation;
      QMouseEvent *mEvent;
      MMoveWire2(doc, 0); // Send the event as 0 so that only drawing update is done

    } else if (QucsMain->MouseMoveAction == &MouseActions::MMoveElement) { // Rotation
      if (selElem && selElem->Type & isComponent) {
        Component *comp = (Component*) selElem;

        if (comp->Ports.count() > 0) { // Do not rotate components without ports
          comp->rotate();
          comp->paintScheme(doc);
        }
      }
    }

    break;

  case Qt::Key_Escape:
    if (QucsMain->MouseMoveAction == &MouseActions::MMoveElement ||
        QucsMain->MouseMoveAction == &MouseActions::MMoveMoving2)
    {
      MReleaseMoving(doc, 0);
    }
    movingElements.clear();
    defaultState(); // Due to "MPressSelect()" capture
    doc->releaseKeyboard();
    //doc->updateContents();
    doc->viewport()->update();
    break;

  /*case Qt::Key_A: //case Qt::Key_Left:
    doc->scrollBy(10, 0);
    break;

  case Qt::Key_D:
    doc->scrollBy(-10, 0);
    break;*/

  default:
    break;
  }

}
//=================================================================================================
/**
 * @brief MouseActions::defaultState Switches to the default "no operation" state.
 * @param Doc
 * @param Event
 */
void MouseActions::defaultState(void) {

  QucsMain->MouseMoveAction = &MouseActions::MMoveFreely;
  QucsMain->MousePressAction = &MouseActions::MPressSelect;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseSelect;
  QucsMain->MouseDoubleClickAction = &MouseActions::MDoubleClickSelect;

  //doc->releaseKeyboard();
  //doc->grabKeyboard(); // Capture all events but do not "accept()" them; sniffing mode
}
//=================================================================================================
/**
 * @brief MouseActions::elementInsertState Switches to the element insertion state.
 * @param Doc
 * @param Event
 */
void MouseActions::elementInsertState(void) {

  QucsMain->MouseMoveAction = &MouseActions::MMoveElement;
  QucsMain->MousePressAction = &MouseActions::MPressElement;
  QucsMain->MouseReleaseAction = &MouseActions::MReleaseElement;
  QucsMain->MouseDoubleClickAction = 0;
}
//=================================================================================================
